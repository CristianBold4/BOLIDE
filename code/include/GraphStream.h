//
// Created by Cristian Boldrin on 24/12/25.
// Reworked for one-pass, bounded-memory, high-throughput sequential parsing.
//

#ifndef GRAPH_STREAM_H
#define GRAPH_STREAM_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>

#if defined(__linux__)
#include <fcntl.h>   // posix_fadvise
#endif

class GraphStream {
public:
    struct Edge {
        size_t u;
        size_t v;
        size_t timestamp;
    };

    // 8 MiB default buffer: good balance for sequential parsing.
    explicit GraphStream(const std::string& filename, size_t buffer_size = 8 * 1024 * 1024)
            : fd_(-1),
              buffer_(buffer_size),
              pos_(0),
              end_(0),
              eof_(false),
              edge_count_(0) {

        if (buffer_size < 1024) {
            throw std::invalid_argument("Buffer size too small");
        }

        fd_ = ::open(filename.c_str(), O_RDONLY);
        if (fd_ == -1) {
            throw std::runtime_error("Failed to open file: " + filename + " - " + std::strerror(errno));
        }

#if defined(__linux__)
        // Hint the kernel that we will read sequentially, once.
        // These are only hints; ignore failures.
        (void) ::posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#ifdef POSIX_FADV_NOREUSE
        (void) ::posix_fadvise(fd_, 0, 0, POSIX_FADV_NOREUSE);
#endif
#endif

        refill();
    }

    ~GraphStream() {
        if (fd_ != -1) {
            ::close(fd_);
        }
    }

    GraphStream(const GraphStream&) = delete;
    GraphStream& operator=(const GraphStream&) = delete;

    size_t getEdgeCount() const noexcept {
        return edge_count_;
    }

    // True if another edge can be parsed.
    // This may pull more data from disk.
    bool hasNext() {
        skipLineBreaksAndSpaces();
        return pos_ < end_ || !eof_;
    }

    bool next(Edge& edge) {
        skipLineBreaksAndSpaces();

        if (pos_ >= end_ && eof_) {
            return false;
        }

        if (!parseUnsigned(edge.u)) {
            return false;
        }

        skipSpaces();

        if (!parseUnsigned(edge.v)) {
            throw std::runtime_error("Malformed input: expected second integer");
        }

        skipSpaces();

        if (!parseUnsigned(edge.timestamp)) {
            throw std::runtime_error("Malformed input: expected third integer");
        }

        skipToNextLine();

        ++edge_count_;
        return true;
    }

private:
    int fd_;
    std::vector<char> buffer_;
    size_t pos_;
    size_t end_;
    bool eof_;
    size_t edge_count_;

    static inline bool isDigit(char c) noexcept {
        return c >= '0' && c <= '9';
    }

    static inline bool isSpaceNoNewline(char c) noexcept {
        return c == ' ' || c == '\t';
    }

    static inline bool isLineBreak(char c) noexcept {
        return c == '\n' || c == '\r';
    }

    void refill() {
        if (eof_) {
            return;
        }

        // Move unread bytes to the beginning.
        if (pos_ > 0 && pos_ < end_) {
            const size_t remaining = end_ - pos_;
            std::memmove(buffer_.data(), buffer_.data() + pos_, remaining);
            end_ = remaining;
            pos_ = 0;
        } else if (pos_ >= end_) {
            pos_ = 0;
            end_ = 0;
        }

        while (end_ < buffer_.size()) {
            const ssize_t n = ::read(fd_, buffer_.data() + end_, buffer_.size() - end_);
            if (n > 0) {
                end_ += static_cast<size_t>(n);
                return;
            }
            if (n == 0) {
                eof_ = true;
                return;
            }
            if (errno == EINTR) {
                continue;
            }
            throw std::runtime_error(std::string("Read failed: ") + std::strerror(errno));
        }

        // Buffer is full and we still need more room.
        // This means one logical token/line is larger than the buffer.
        // That should never happen for normal edge-list files.
        if (pos_ == 0 && end_ == buffer_.size()) {
            throw std::runtime_error("Buffer too small: encountered token/line larger than buffer");
        }
    }

    inline void ensureData() {
        if (pos_ >= end_ && !eof_) {
            refill();
        }
    }

    void skipSpaces() {
        while (true) {
            while (pos_ < end_ && isSpaceNoNewline(buffer_[pos_])) {
                ++pos_;
            }
            if (pos_ < end_ || eof_) {
                return;
            }
            refill();
        }
    }

    void skipLineBreaksAndSpaces() {
        while (true) {
            while (pos_ < end_) {
                const char c = buffer_[pos_];
                if (isSpaceNoNewline(c) || isLineBreak(c)) {
                    ++pos_;
                } else {
                    return;
                }
            }
            if (eof_) {
                return;
            }
            refill();
        }
    }

    void skipToNextLine() {
        while (true) {
            while (pos_ < end_) {
                const char c = buffer_[pos_++];
                if (c == '\n') {
                    return;
                }
            }
            if (eof_) {
                return;
            }
            refill();
        }
    }

    bool parseUnsigned(size_t& value) {
        value = 0;
        bool saw_digit = false;

        while (true) {
            while (pos_ < end_) {
                const char c = buffer_[pos_];
                if (!isDigit(c)) {
                    return saw_digit;
                }
                saw_digit = true;
                value = value * 10 + static_cast<size_t>(c - '0');
                ++pos_;
            }

            if (eof_) {
                return saw_digit;
            }

            refill();
        }
    }
};

#endif // GRAPH_STREAM_H