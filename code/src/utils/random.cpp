//
// Created by Cristian Boldrin on 14/11/25.
//

#include <utils/random.h>

using namespace utils;

// static std::mt19937 randomEngine;

RandomIndexer::RandomIndexer(size_t s, std::mt19937 &engine) : sampler(0, s - 1), idxEngine_(engine)
{
}

size_t
RandomIndexer::next()
{
    return sampler(idxEngine_);
}

void
Random::initialize(size_t fixedSeed)
{
    if (fixedSeed == -1)
    {
        // Source: https://stackoverflow.com/questions/15509270/does-stdmt19937-require-warmup
        // -- warmup std mt
        std::vector<std::random_device::result_type> seedData( std::mt19937::state_size );
        std::random_device randomDevice;
        std::generate_n(seedData.data(), seedData.size(), std::ref(randomDevice));
        std::seed_seq randomSeq(std::begin(seedData), std::end(seedData));
        randomEngine_.seed(randomSeq);
    }
    else
    {
        randomEngine_.seed(static_cast<uint>(fixedSeed));
    }
}

RandomIndexer
Random::getIndexer(size_t size)
{
    return RandomIndexer(size, randomEngine_);
}

double
Random::getDouble()
{
    return pickRandomValue(randomEngine_);
    // return randomEngine_() * (1.0 / 4294967296.0); // faster but only 32 bits of randomness
}

Random::Random()
{
}

Random::Random(size_t fixedSeed)
{
    Random::initialize(fixedSeed);
}

size_t
Random::stochasticRounding(double value)
{
    auto valueHigh = floor(value);
    auto valueLow = ceil(value);
    auto proba = (value - valueLow) / (valueHigh - valueLow);
    auto randomVal = this->getDouble();
    if (randomVal < proba)
    {
        return static_cast<size_t>(round(valueHigh)); // Round up
    }
    return static_cast<size_t>(round(valueLow)); // Round down
}
