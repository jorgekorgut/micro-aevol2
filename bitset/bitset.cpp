#include "bitset.h"
#include <cstring>
#include <algorithm>
#include <iostream>
#include <cassert>

Bitset::Bitset(size_t bitsetSize)
{
    blockSizeBytes = sizeof(u_int64_t);
    blockSizeBites = blockSizeBytes * 8;
    size = bitsetSize;
    blockCount = (size / (blockSizeBites) + 1) + 2;
    blocks = new u_int64_t[blockCount];

    std::memset(blocks, 0, sizeof(u_int64_t) * (blockCount));
}

Bitset::Bitset(const Bitset &clone)
{
    blockSizeBytes = sizeof(u_int64_t);
    blockSizeBites = blockSizeBytes * 8;
    size = clone.bitsetSize();
    blockCount = (size / (blockSizeBites) + 1) + 2;
    blocks = new u_int64_t[blockCount];

    std::memcpy(blocks, clone.getBlocks(), blockSizeBytes * blockCount);
}

Bitset::Bitset(const std::string &bitset, size_t bitsetSize, bool reverseString)
{
    blockSizeBytes = sizeof(u_int64_t);
    blockSizeBites = blockSizeBytes * 8;
    size = bitsetSize;
    blockCount = (size / (blockSizeBites) + 1) + 2;
    blocks = new u_int64_t[blockCount];

    std::memset(blocks, 0, blockSizeBytes * blockCount);

    std::string currentBitset = bitset;

    u_int64_t iteratorMask = 1;

    if (reverseString)
    {
        std::reverse(currentBitset.begin(), currentBitset.end());
    }

    size_t globalBitIterator = 0;
    int currentBlockIterator = 1;
    int currentIterator = 0;
    u_int64_t *currentBlock = blocks + currentBlockIterator;

    while (globalBitIterator < size)
    {
        if (currentBitset[globalBitIterator] == '1')
        {
            *currentBlock |= iteratorMask;
        }
        iteratorMask <<= 1;

        ++currentIterator;
        ++globalBitIterator;
        if (currentIterator >= blockSizeBites)
        {
            iteratorMask = 1;
            currentIterator = 0;
            ++currentBlockIterator;
            currentBlock = blocks + currentBlockIterator;
        }
    }
}

std::string Bitset::printDebug() const
{
    return printDebug(blocks, blockSizeBites);
}

std::string Bitset::printDebug(u_int64_t *localBlocks, int localBlockSizeBites) const
{
    std::string response = "";
    u_int64_t iteratorMask = 1;

    for (size_t currentBlockIterator = 0; currentBlockIterator < blockCount; ++currentBlockIterator)
    {
        u_int64_t currentElement = localBlocks[currentBlockIterator];
        for (size_t j = 0; j < localBlockSizeBites; ++j)
        {
            response += ('0' + (currentElement & iteratorMask));
            currentElement >>= 1;
        }
    }

    std::cout << response;
    std::cout << " | block offset : " << insideBlockOffset;
    std::cout << std::endl;
    return response;
}

std::string Bitset::print()
{
    std::string response = "";
    u_int64_t iteratorMask = 1;
    int currentBlockIterator = 1;
    int currentIterator = insideBlockOffset;
    size_t globalBitIterator = 0;
    if (insideBlockOffset < 0)
    {
        currentBlockIterator = 0;
        currentIterator = blockSizeBites + insideBlockOffset;
    }

    u_int64_t currentBlock = blocks[currentBlockIterator];

    if (insideBlockOffset > 0)
    {
        currentBlock >>= insideBlockOffset;
    }

    while (globalBitIterator < size)
    {
        response += ('0' + (currentBlock & iteratorMask));
        currentBlock >>= 1;

        ++currentIterator;
        ++globalBitIterator;
        if (currentIterator >= blockSizeBites)
        {
            currentIterator = 0;
            ++currentBlockIterator;
            currentBlock = blocks[currentBlockIterator];
        }
    }
    std::cout << response << std::endl;
    return response;
}

Bitset &Bitset::operator>>=(int shiftNumber)
{
    int offsetBlockShift = (insideBlockOffset - shiftNumber) / blockSizeBites; // Take care of offsetBlockShift int overflow
    insideBlockOffset = (insideBlockOffset - shiftNumber) % blockSizeBites;

    std::memcpy(blocks, blocks + offsetBlockShift, (blockCount + offsetBlockShift) * blockSizeBytes);

    std::memset(blocks, 0, (-offsetBlockShift) * blockSizeBytes);
    /*
    for (int currentBlock = blockCount - 1; currentBlock + offsetBlockShift >= 0; currentBlock--)
    {
        // std::cout << currentBlock << "<-" << currentBlock + offsetBlockShift << std::endl;
        blocks[currentBlock] = blocks[currentBlock + offsetBlockShift];
    }

    for (size_t currentBlock = 0; currentBlock < -offsetBlockShift; currentBlock++)
    {
        blocks[currentBlock] = 0;
    }
    */
    return *this;
}

Bitset &Bitset::operator<<=(int shiftNumber)
{
    int offsetBlockShift = (insideBlockOffset + shiftNumber) / blockSizeBites; // Take care of offsetBlockShift int overflow
    insideBlockOffset = (insideBlockOffset + shiftNumber) % blockSizeBites;

    std::memcpy(blocks, blocks + offsetBlockShift, (blockCount - offsetBlockShift) * blockSizeBytes);

    std::memset(blocks - (blockCount - offsetBlockShift), 0, (offsetBlockShift)*blockSizeBytes);
    /*
    for (int currentBlock = blockCount - 1; currentBlock + offsetBlockShift >= 0; currentBlock--)
    {
        // std::cout << currentBlock << "<-" << currentBlock + offsetBlockShift << std::endl;
        blocks[currentBlock] = blocks[currentBlock + offsetBlockShift];
    }

    for (size_t currentBlock = 0; currentBlock < -offsetBlockShift; currentBlock++)
    {
        blocks[currentBlock] = 0;
    }
    */
    return *this;
}

bool Bitset::compare(int fromIndex, const Bitset &compareTo, int toIndex, int length) const
{
    int lengthBlock = length / blockSizeBites + 1;
    u_int64_t *copyCompareFromBlocks = new u_int64_t[lengthBlock];
    u_int64_t *copyCompareToBlocks = new u_int64_t[lengthBlock];

    alignBlocks(copyCompareFromBlocks, fromIndex, copyCompareToBlocks, compareTo, toIndex, length, lengthBlock);

    doReset(copyCompareFromBlocks, blockSizeBites, length, lengthBlock * blockSizeBites);
    doReset(copyCompareToBlocks, compareTo.getBlockSizeBites(), length, lengthBlock * compareTo.getBlockSizeBites());

    // std::cout << "Reseted blocks from bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // compareTo.printDebug(copyCompareToBlocks, blockSizeBites);

    int response = memcmp(copyCompareToBlocks, copyCompareFromBlocks, lengthBlock * blockSizeBytes);

    delete[] copyCompareFromBlocks;
    delete[] copyCompareToBlocks;

    return (response == 0);

    // int globalIterator = 0;
    // u_int64_t mask = 1;

    // int blockFromIndex = fromIndex / blockSizeBites + 1;
    // int bitFromIndex = fast_mod(fromIndex, blockSizeBites);
    // u_int64_t currentFromBlock = blocks[blockFromIndex];
    // currentFromBlock >>= bitFromIndex;

    // int blockToIndex = toIndex / blockSizeBites + 1;
    // int bitToIndex = fast_mod(toIndex, blockSizeBites);
    // u_int64_t currentToBlock = compareTo.getBlocks()[blockToIndex];
    // currentToBlock >>= bitToIndex;

    // int isEqual = true;
    // while (isEqual && globalIterator < length)
    // {
    //     if((mask & currentFromBlock) == (mask & currentToBlock)){
    //         isEqual = false;
    //         break;
    //     }
    //     currentFromBlock >>= 1;
    //     currentToBlock >>= 1;

    //     ++bitToIndex;
    //     if (bitToIndex >= blockSizeBites)
    //     {
    //         ++blockToIndex;
    //         currentToBlock = compareTo.getBlocks()[blockToIndex];
    //         bitToIndex = 0;
    //     }

    //     ++bitFromIndex;
    //     if (bitFromIndex >= blockSizeBites)
    //     {
    //         ++blockFromIndex;
    //         currentFromBlock = blocks[blockFromIndex];
    //         bitFromIndex = 0;
    //     }

    //     ++globalIterator;
    // }
    // return isEqual;
}

bool Bitset::compareIgnore(int fromIndex, const Bitset &compareTo, int toIndex, int length, int ignoreIndex, int ignoreLength) const
{
    int lengthBlock = length / blockSizeBites + 1;
    u_int64_t *copyCompareFromBlocks = new u_int64_t[lengthBlock];
    u_int64_t *copyCompareToBlocks = new u_int64_t[lengthBlock];

    alignBlocks(copyCompareFromBlocks, fromIndex, copyCompareToBlocks, compareTo, toIndex, length, lengthBlock);

    doReset(copyCompareFromBlocks, blockSizeBites, length, lengthBlock * blockSizeBites);
    // doReset(copyCompareFromBlocks, blockSizeBites, ignoreIndex, ignoreIndex + ignoreLength - 1);
    doReset(copyCompareToBlocks, compareTo.getBlockSizeBites(), length, lengthBlock * compareTo.getBlockSizeBites());
    // doReset(copyCompareToBlocks, compareTo.getBlockSizeBites(), ignoreIndex, ignoreIndex + ignoreLength - 1);

    // std::cout << "Reseted blocks from bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // compareTo.printDebug(copyCompareToBlocks, blockSizeBites);

    int response = memcmp(copyCompareToBlocks, copyCompareFromBlocks, lengthBlock * blockSizeBytes);

    delete[] copyCompareFromBlocks;
    delete[] copyCompareToBlocks;

    return (response == 0);
}

bool Bitset::containsPatternAtPosition(int position, const Bitset &pattern, int patternLength) const
{
    bool response = true;
    int blockIndex = position / blockSizeBites + 1;
    int bitIndex = fast_mod(position, blockSizeBites);

    if (bitIndex + pattern.bitsetSize() >= blockSizeBites)
    {
        // u_int64_t leftMask = ~0;
        // leftMask <<= bitIndex;
        // leftMask = ~leftMask;
        // u_int64_t patternMaskLeft = leftMask & (pattern.getBlocks()[0] >> bitIndex);

        // u_int64_t rightMask = ~0;
        // rightMask >>= pattern.getBlockSizeBites() - bitIndex - 1;
        // rightMask = ~rightMask;
        // u_int64_t patternMaskRight = rightMask & (pattern.getBlocks()[0] >> pattern.getBlockSizeBites() - bitIndex - 1);

        // response &= (leftMask & blocks[blockIndex] == patternMaskLeft);
        // response &= (rightMask & blocks[blockIndex + 1] == patternMaskRight);

        int t_pos = position;
        for (int k = 0; k < patternLength; k++, t_pos++)
        {
            if (t_pos >= size)
            {
                t_pos -= size;
            }
            if (operator[](t_pos) != pattern[k])
            {
                response = false;
                break;
            }
        }
    }
    else
    {

        u_int64_t leftMask = ~0;
        leftMask <<= bitIndex;
        u_int64_t rightMask = ~0;
        rightMask <<= bitIndex + patternLength;
        rightMask = ~rightMask;

        u_int64_t patternMaskLeft = pattern.getBlocks()[1] << bitIndex;

        response = ((leftMask & blocks[blockIndex] & rightMask) == patternMaskLeft);
    }
    return response;
}

bool Bitset::containsPatternAtPositionIgnore(int position, const Bitset &pattern, int patternLength, int ignorePos, int ignoreLength, int postPatternSize) const
{
    bool response = true;
    int blockIndex = position / blockSizeBites + 1;
    int bitIndex = fast_mod(position, blockSizeBites);

    if (bitIndex + patternLength + ignoreLength + postPatternSize >= blockSizeBites)
    {

        u_int64_t leftMask = ~0;
        leftMask <<= bitIndex;
        u_int64_t patternMaskLeft = leftMask & (pattern.getBlocks()[1] << bitIndex);

        response = (leftMask & blocks[blockIndex] == patternMaskLeft);
        if (response == false)
        {
            return response;
        }

        leftMask = ~0;
        leftMask <<= blockSizeBites - (bitIndex + patternLength);
        leftMask = ~leftMask;

        u_int64_t patternMaskRight = leftMask & (pattern.getBlocks()[1] >> bitIndex);

        response = (leftMask & blocks[blockIndex] == patternMaskRight);

        if (response == false)
        {
            return response;
        }

        int t_pos = position + patternLength + ignoreLength;
        int bitsetSize = size;

        for (int k = 0; k < postPatternSize; ++k, ++t_pos)
        {
            if (t_pos >= bitsetSize)
                t_pos -= bitsetSize;

            if (operator[](t_pos) != pattern[k + patternLength + ignoreLength])
            {
                response = false;
                break;
            }
        }
    }
    else
    {
        u_int64_t leftMask = ~0;
        leftMask <<= bitIndex;
        u_int64_t rightMask = ~0;
        rightMask <<= bitIndex + patternLength;
        rightMask = ~rightMask;

        u_int64_t patternMaskLeft = pattern.getBlocks()[1] << bitIndex;
        response = ((leftMask & blocks[blockIndex] & rightMask) == (patternMaskLeft));
        if (response == false)
        {
            return response;
        }

        rightMask = ~0;
        rightMask <<= bitIndex + patternLength + ignoreLength + postPatternSize;
        rightMask = ~rightMask;

        u_int64_t ignoreLeftMask = ~0;
        ignoreLeftMask <<= bitIndex + patternLength;
        u_int64_t ignoreRightMask = ~0;
        ignoreRightMask <<= bitIndex + patternLength + ignoreLength;
        ignoreRightMask = ~ignoreRightMask;

        u_int64_t ignoreMask = ~(ignoreLeftMask & ignoreRightMask);
        response = ((ignoreMask & leftMask & blocks[blockIndex] & rightMask) == (ignoreMask & patternMaskLeft));
    }
    return response;
}

// FIXME: Negative numbers are not taken into account
// u_int64_t Bitset::getMask(int pos, int length)
// {
//     int correctedIndex = fast_mod(pos, size);
//     u_int64_t value = 0;
//     int blockIndex = correctedIndex / blockSizeBites + 1;

//     if (blockIndex > blockCount - 2)
//     {
//         blockIndex = 1;
//     }

//     int nextBlockIndex = blockIndex + 1;

//     if (nextBlockIndex > blockCount - 2)
//     {
//         nextBlockIndex = 1;
//     }

//     int bitIndex = fast_mod(correctedIndex, blockSizeBites);

//     // If the mask is truncated by blocks or end of bitset
//     if (bitIndex + length >= blockSizeBites || pos + length >= size)
//     {
//         u_int64_t leftMask = ~0;
//         leftMask <<= length;
//         leftMask = ~leftMask;

//         value = blocks[blockIndex] >> bitIndex | blocks[nextBlockIndex] << (blockSizeBites - bitIndex);

//         if (pos + length >= size)
//         {
//             int lastBlockRealSize = fast_mod(size, blockSizeBites);
//             value |= blocks[nextBlockIndex] << (blockSizeBites - bitIndex) + lastBlockRealSize;
//         }

//         value = value & leftMask;

//     }
//     else
//     {
//         u_int64_t leftMask = ~0;
//         leftMask <<= bitIndex;
//         u_int64_t rightMask = ~0;
//         rightMask <<= bitIndex + length;
//         rightMask = ~rightMask;

//         value = (leftMask & blocks[blockIndex] & rightMask) >> bitIndex;
//     }

//     return value;
// }

int Bitset::compareDistance(int fromIndex, const Bitset &compareTo, int toIndex, int length) const
{
    int lengthBlock = length / blockSizeBites + 1;
    u_int64_t *copyCompareFromBlocks = new u_int64_t[lengthBlock];
    u_int64_t *copyCompareToBlocks = new u_int64_t[lengthBlock];

    alignBlocks(copyCompareFromBlocks, fromIndex, copyCompareToBlocks, compareTo, toIndex, length, lengthBlock);

    // std::cout << "Aligned blocks from bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // compareTo.printDebug(copyCompareToBlocks, blockSizeBites);

    int currentBlockIterator = 0;
    int globalBitIterator = 0;
    int currentIterator = 0;

    u_int64_t currentBlockFrom = copyCompareFromBlocks[currentBlockIterator];
    u_int64_t currentBlockTo = copyCompareToBlocks[currentBlockIterator];

    u_int64_t iteratorMask = 1;
    int count = 0;

    while (globalBitIterator < length)
    {
        if ((currentBlockFrom & iteratorMask) != (currentBlockTo & iteratorMask))
        {
            ++count;
        }
        iteratorMask <<= 1;
        ++globalBitIterator;
        ++currentIterator;
        if (currentIterator >= blockSizeBites)
        {
            iteratorMask = 1;
            currentIterator = 0;
            ++currentBlockIterator;
            currentBlockFrom = copyCompareFromBlocks[currentBlockIterator];
            currentBlockTo = copyCompareToBlocks[currentBlockIterator];
        }
    }

    delete[] copyCompareFromBlocks;
    delete[] copyCompareToBlocks;

    return count;
}

/*
 * fromIndex in virtual space
 */
void Bitset::alignBlocks(u_int64_t *copyCompareFromBlocks, int fromIndex, u_int64_t *copyCompareToBlocks, const Bitset &compareTo, int toIndex, int length, int lengthBlock) const
{
    u_int64_t *const compareFromBlocks = blocks;
    u_int64_t *const compareToBlocks = compareTo.getBlocks();

    /*
    std::cout << "Enter bitsets: " << std::endl;
    this->printDebug(compareFromBlocks, blockSizeBites);
    compareTo.printDebug(compareToBlocks, compareTo.getBlockSizeBites());
    */

    toIndex += compareTo.getInsideBlockOffset();
    int toBlockIndex = toIndex / compareTo.getBlockSizeBites() + 1;
    int toBitIndex = fast_mod(toIndex, compareTo.getBlockSizeBites());
    int toLengthBlock = lengthBlock;

    if (toBitIndex < 0)
    {
        toBlockIndex -= 1;
        toBitIndex = compareTo.getBlockSizeBites() + toBitIndex;
        ++toLengthBlock;
    }

    fromIndex += insideBlockOffset;
    int fromBlockIndex = fromIndex / blockSizeBites + 1;
    int fromBitIndex = fast_mod(fromIndex, blockSizeBites);
    int fromLengthBlock = lengthBlock;

    if (fromBitIndex < 0)
    {
        fromBlockIndex -= 1;
        fromBitIndex = blockSizeBites + fromBitIndex;
        ++fromLengthBlock;
    }

    // std::cout << "fromBlockIndex:" << fromBlockIndex << std::endl;
    // std::cout << "fromBitIndex:" << fromBitIndex << std::endl;
    // std::cout << "FromInsideBlockOffset:" << insideBlockOffset << std::endl;
    // std::cout << "fromLengthBlock:" << fromLengthBlock << std::endl;
    // std::cout << "toBlockIndex:" << toBlockIndex << std::endl;
    // std::cout << "toBitIndex:" << toBitIndex << std::endl;
    // std::cout << "lengthBlock:" << lengthBlock << std::endl;

    memcpy(copyCompareFromBlocks, compareFromBlocks + fromBlockIndex, fromLengthBlock * blockSizeBytes);

    memcpy(copyCompareToBlocks, compareToBlocks + toBlockIndex, toLengthBlock * compareTo.getBlockSizeBytes());

    // std::cout << "Copied bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // this->printDebug(copyCompareToBlocks, blockSizeBites);

    // Shift few elements (boost/dynamic_bitset)
    if (fromBitIndex > 0)
    {
        const int leftBitShift = blockSizeBites - fromBitIndex;

        for (int i = 0; i < fromLengthBlock - 1; ++i)
        {
            copyCompareFromBlocks[i] = (copyCompareFromBlocks[i] >> fromBitIndex) | (copyCompareFromBlocks[i + 1] << leftBitShift);
        }

        copyCompareFromBlocks[fromLengthBlock - 1] = copyCompareFromBlocks[fromLengthBlock - 1] >> fromBitIndex;
    }

    if (toBitIndex > 0)
    {
        const int leftBitShift = compareTo.getBlockSizeBites() - toBitIndex;

        for (int i = 0; i < toLengthBlock - 1; ++i)
        {
            copyCompareToBlocks[i] = (copyCompareToBlocks[i] >> toBitIndex) | (copyCompareToBlocks[i + 1] << leftBitShift);
        }

        copyCompareToBlocks[toLengthBlock - 1] = copyCompareToBlocks[toLengthBlock - 1] >> toBitIndex;
    }

    // std::cout << "Copied blocks from bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // compareTo.printDebug(copyCompareToBlocks, blockSizeBites);
}

// fromBitIndex inclusive, ToBitIndex exclusive
// It doesnt take into account right and left offset inside block
void Bitset::doReset(u_int64_t *blocks, int blockBiteSize, int fromBiteIndex, int toBiteIndex)
{
    u_int64_t constMask = 1;

    int distance = toBiteIndex - fromBiteIndex;
    int completeBlocksSize = distance / blockBiteSize;

    int currentBlockBiteIndex = fast_mod(fromBiteIndex, blockBiteSize);
    int currentBlockIndex = fromBiteIndex / blockBiteSize;

    u_int64_t *currentBlock = blocks + currentBlockIndex;

    // std::cout << "distance:" << distance << std::endl;
    // std::cout << "completeBlocksSize:" << completeBlocksSize << std::endl;
    // std::cout << "currentBlockIndex:" << currentBlockIndex << std::endl;
    // std::cout << "currentBlockBiteIndex:" << currentBlockBiteIndex << std::endl;

    constMask <<= currentBlockBiteIndex;

    while (fromBiteIndex < toBiteIndex)
    {

        (*currentBlock) &= (~constMask);

        constMask <<= 1;
        ++fromBiteIndex;
        ++currentBlockBiteIndex;
        if (currentBlockBiteIndex >= blockBiteSize)
        {
            memset(currentBlock, 0, completeBlocksSize * (blockBiteSize / 8));
            currentBlock += completeBlocksSize;

            currentBlockBiteIndex = 0;
            constMask = 1;
        }
    }
}

bool Bitset::operator[](int targetNumber) const
{
    targetNumber += insideBlockOffset;
    int currentBlockIndex = targetNumber / blockSizeBites + 1;
    int currentBlockBiteIndex = fast_mod(targetNumber, blockSizeBites);

    return (blocks[currentBlockIndex] >> currentBlockBiteIndex) & 1;
}

void Bitset::set(int targetNumber, bool value)
{
    targetNumber += insideBlockOffset;
    int currentBlockIndex = targetNumber / blockSizeBites + 1;
    int currentBlockBiteIndex = fast_mod(targetNumber, blockSizeBites);
    u_int64_t removeMask = 1;
    removeMask <<= currentBlockBiteIndex;

    if (value == false)
    {
        blocks[currentBlockIndex] &= (~removeMask);
    }
    else
    {
        blocks[currentBlockIndex] |= (removeMask);
    }
}

void Bitset::flip(int targetNumber)
{
    targetNumber += insideBlockOffset;
    int currentBlockIndex = targetNumber / blockSizeBites + 1;
    int currentBlockBiteIndex = fast_mod(targetNumber, blockSizeBites);
    u_int64_t flipMask = 1;

    blocks[currentBlockIndex] ^= flipMask << currentBlockBiteIndex;
    /*
    for (int32_t i = 0; i < bitsetSize(); i++)
    {
        std::cerr << this->operator[](i);
    }
    std::cerr << std::endl;
    */
}

int Bitset::getInsideBlockOffset() const
{
    return insideBlockOffset;
}

size_t Bitset::bitsetSize() const
{
    return size;
}

u_int64_t *Bitset::getBlocks() const
{
    return blocks;
}

int Bitset::getBlockSizeBites() const
{
    return blockSizeBites;
}

int Bitset::getBlockSizeBytes() const
{
    return blockSizeBytes;
}

int Bitset::getBlockCount() const
{
    return blockCount;
}

Bitset::~Bitset()
{
    // delete[] blocks;
}