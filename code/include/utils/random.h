//
// Created by Cristian Boldrin on 14/11/25.
//

#pragma once

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <unordered_set>

namespace utils {

    class RandomIndexer {
    public:
        RandomIndexer(size_t size, std::mt19937 &engine);

        size_t next();

    private:
        std::uniform_int_distribution<size_t> sampler;
        std::mt19937& idxEngine_;
    };

    class Random {
    public:
        /**
         * @brief Initialises random class.
         */
        Random();

        /**
         * @brief Initialises random class with a fixed seed.
         * @param fixedSeed
         */
        Random(size_t fixedSeed);

        void initialize(size_t fixedSeed);

        /**
         * Returns a random real number in the interval [0.0, 1.0).
         */
        double
        getDouble();

        /**
        * Returns a random real number in the interval [min, max).
        */
        double
        getDouble(double min, double max);

        std::mt19937 & getEngine() { return randomEngine_; }

        RandomIndexer
        getIndexer(size_t size);


        /**
         * @brief Select a number of elements from vector uniformly at random.
         * @param elements A collection of elements to sample from.
         * @param numberOfElements The number of samples to return.
         */
        template<typename T>
        std::vector<T>
        choice(const std::vector<T> &elements, const size_t numberOfElements) {
            // Notice that the templated class function is implemented in this header
            // file because if we tried to move it in the source file then the linker
            // error with "undefined reference to". For more information about this issue
            // read https://bytefreaks.net/programming-2/c/c-undefined-reference-to-templated-class-function
            std::vector<T> samples;
            auto indexSampler = getIndexer(elements.size());

            for (size_t i = 0; i < numberOfElements; i++) {
                auto sampledIndex = indexSampler.next();
                samples.push_back(elements[sampledIndex]);
            }

            return samples;
        }

        /**
         * @brief Stochastically rounds up or down a real number `v` with probability (v-⌊v⌋)/(⌈v⌉-⌊v⌋).
         * See https://nhigham.com/2020/07/07/what-is-stochastic-rounding/
         * @param value The floating point number to be round up or down.
         */
        size_t
        stochasticRounding(double value);


    private:

        std::mt19937 randomEngine_;
        std::uniform_real_distribution<> pickRandomValue;

    };

}

