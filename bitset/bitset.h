#pragma once

#include <string>

inline int fast_mod(const int input, const int ceil)
{
    // apply the modulo operator only when needed
    // (i.e. when the input is greater than the ceiling)
    return input >= ceil ? input % ceil : input;
    // NB: the assumption here is that the numbers are positive
}

class Bitset
{
public:
    Bitset() = default;
    Bitset(size_t size);
    Bitset(const std::string &bitset, size_t size, bool reverseString = false);
    Bitset(const Bitset &clone);
    ~Bitset();
    std::string printDebug(u_int64_t *localBlocks, int localBlockSize) const;
    std::string printDebug() const;
    std::string print();
    size_t bitsetSize() const;
    void set(int targetNumber, bool value);
    void flip(int targetNumber);
    u_int64_t getMask(int pos, int length);
    bool compare(int fromIndex, const Bitset &compareTo, int toIndex, int length) const;
    bool compareIgnore(int fromIndex, const Bitset &compareTo, int toIndex, int length, int ignoreIndex, int ignoreLength) const;
    bool containsPatternAtPosition(int position, const Bitset &pattern, int patternLength) const;
    bool containsPatternAtPositionIgnore(int position, const Bitset &pattern, int patternLength, int ignorePos, int ignoreLength, int postPatternSize) const;
    int compareDistance(int fromIndex, const Bitset &compareTo, int toIndex, int length) const;

    Bitset &operator>>=(int shiftNumber);
    Bitset &operator<<=(int shiftNumber);
    bool operator[](int targetNumber) const;

    u_int64_t *getBlocks() const;
    int getBlockSizeBites() const;
    int getBlockSizeBytes() const;
    int getBlockCount() const;
    static void doReset(u_int64_t *blocks, int blockBiteSize, int from, int to);
    int getInsideBlockOffset() const;
    void alignBlocks(u_int64_t *copyCompareFromBlocks, int fromIndex, u_int64_t *copyCompareToBlocks, const Bitset &compareTo, int toIndex, int length, int lengthBlock) const;

private:
    u_int64_t *blocks = nullptr;
    // Number of bits in the bitset
    int size = 0;

    // number of bytes in 1 block
    const static int blockSizeBytes = sizeof(u_int64_t);

    // number of bites in 1 block (blockSizeBytes * 8)
    const static int blockSizeBites = blockSizeBytes * 8;

    // Number of blocks in the bitset
    int blockCount = 0;

    int insideBlockOffset = 0;
};