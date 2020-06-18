/*
 *
 *  Public Domain Zip - Public domain ZIP archiver implementation
 *  Written in 2018 by Martin Rizzo <martinrizzo@gmail.com>
 *
 *  To the extent possible under law, the author(s) have dedicated all
 *  copyright and related and neighboring rights to this software to
 *  the public domain worldwide. This software is distributed WITHOUT
 *  ANY WARRANTY.
 *
 *  You should have received a copy of the CC0 Public Domain Dedication
 *  along with this software. If not, see:
 *  http://creativecommons.org/publicdomain/zero/1.0/
 *
 */
#include <cstdio>
#include <cstring>
#include "PDZip_ReversedHuffmanDecoder.h"

//======================================================================================================================
#   pragma mark - INTERNAL LINKAGE [only for this transaction unit]

#define Data        PDZip::ReversedHuffmanDecoder::Data
#define CodeLengths PDZip::CodeLengths


static inline unsigned reversedINC(unsigned huffman, unsigned length) {
    assert( length>0 );
    unsigned incBit = 1<<(length-1);
    do { huffman ^= incBit; } while ( (huffman&incBit)==0 && (incBit>>=1)!=0 );
    return huffman;
}

static inline void setWithRepetitions(Data*    table,
                                      unsigned availableEntries,
                                      unsigned huffman,
                                      Data     data,
                                      unsigned length)
{
    assert( table!=NULL );
    assert( huffman < (1<<length) );
    assert( length>=1 );
    const unsigned unknownStep = (1<<length);
    for (unsigned unknownBits=0; unknownBits<availableEntries; unknownBits+=unknownStep) {
        table[ unknownBits | huffman ] = data;
    }
}

static inline int makeSubTable(Data*              subtable,
                               int                availableEntries,
                               unsigned           firstHuffman,
                               unsigned           maxLength,
                               const CodeLengths& codeLengths,
                               int                firstIndex,
                               int                endIndex)
{
    static const unsigned DiscardedBits = 8; // < number of bits discarded by the subtable
    const int numberOfEntries = 1<<(maxLength-DiscardedBits);
    
    assert( numberOfEntries<=availableEntries  );
    
    unsigned huffman = firstHuffman;
    int      index   = firstIndex;
    while ( index!=endIndex ) {
        const unsigned code   = codeLengths.getCodeAt(index);
        const unsigned length = codeLengths.getLengthAt(index);
        const Data     data   = Data::fromCodeLength(code,length);
        setWithRepetitions(subtable, numberOfEntries, huffman>>DiscardedBits, data, length-DiscardedBits);

        huffman = reversedINC(huffman,length);
        index   = codeLengths.getNextIndex(index);
    }
    return numberOfEntries;
}

#undef Data
#undef CodeLengths


namespace PDZip {

    
    //==================================================================================================================
#   pragma mark - LOADING THE HUFFMAN TABLE

    
    void ReversedHuffmanDecoder::load(const CodeLengths& codeLengths) {
        assert( codeLengths.isOpen()==false );
        
        // reset the main-table
        std::memset(_table, 0, MainTableSize*sizeof(Data));
        
        unsigned huffman = 0;
        int      index   = codeLengths.getFirstIndex();
        unsigned length  = codeLengths.getLengthAt(index);

        // lengths from 1 to 8
        // unknown bits are filled with all possible values
        while ( index!=0 && length<=8 ) {
            const Data data = Data::fromCodeLength(codeLengths.getCodeAt(index),length);
            setWithRepetitions(_table, 256, huffman, data, length);
            huffman = reversedINC(huffman,length);
            index   = codeLengths.getNextIndex(index);
            length  = codeLengths.getLengthAt(index);
        }
        // lengths from 9
        // subtables are created
        int insertIndex = MainTableSize;
        while ( index!=0 ) {
            const unsigned firstHuffman = huffman;
            const int      firstIndex   = index;
            do {
                length  = codeLengths.getLengthAt(index);
                huffman = reversedINC(huffman,length);
                index   = codeLengths.getNextIndex(index);
            } while ( index!=0 && (huffman&0xFF)==(firstHuffman&0xFF) );
            
            const unsigned maxLength = length;
            _table[firstHuffman] = Data::fromIndexMask(insertIndex, (1<<(maxLength-8))-1);
            insertIndex += makeSubTable(&_table[insertIndex],
                                        (FullTableSize-insertIndex),
                                        firstHuffman, maxLength,
                                        codeLengths, firstIndex, index
                                        );
        }
        //static bool tableWasPrinted = false; if (!tableWasPrinted) { tableWasPrinted = true;
        //    _printDebugInfo();
        //}
    }

    
    //==================================================================================================================
#   pragma mark - GETTING DEBUG INFORMATION

#   if !defined(NDEBUG)

    /**
     * Prints to console debug information about the current state of this object
     */
    void CodeLengths::_printDebugInfo() {

        std::printf("(SortedCodeLengths) {\n");
        if (isOpen()) {
            std::printf(" ... is still open ... ");
        }
        else {
            unsigned lastLength       = 0xFFFFFF;
            for (int index = getFirstIndex(); index!=0; index = getNextIndex(index)) {
                unsigned length = getLengthAt(index);
                if (lastLength!= length) {
                    if (lastLength!=0xFFFFFF) { std::printf("}\n"); }
                    std::printf(" Length %d: {", length);
                    lastLength = length;
                }
                else {
                    std::printf(", ");
                }
                std::printf("%d", getCodeAt(index));
            }
        }
        std::printf("}\n");
    }
    
    /**
     * Prints to console the current state of the internal huffman table
     */
    void ReversedHuffmanDecoder::printDebugInfo() {
        bool isPrintedFlags[256];
        for (int i=0; i<256; ++i) { isPrintedFlags[i]=false; }
        
        std::printf("\n");
        unsigned index=0; while (index<256) {
            if ( isPrintedFlags[index]==false ) {
                if      ( _table[index].raw()==0  ) { _printDebugInfo_TableEntry(_table[index], index, 0); isPrintedFlags[index]=true; }
                else if ( _table[index].isValid() ) { _printDebugInfo_CodeRepetition(_table, index, isPrintedFlags); }
                else                                { _printDebugInfo_Subtable(_table, index, isPrintedFlags); }
            }
            index = index<255 ? reversedINC(index,8) : 256;
        }
        std::printf("\n");
    }

    /**
     * Prints debugging information about "code,length" repetitions in the internal huffman table
     */
    void ReversedHuffmanDecoder::_printDebugInfo_CodeRepetition(Data* table, int firstIndex, bool* isPrintedFlags) {
        assert( table!=NULL && firstIndex<MainTableSize && isPrintedFlags!=NULL );
        assert( table[firstIndex].isValid() );
        const Data     data             = table[firstIndex];
        const unsigned rawRepeatedValue = data.raw();
        
        // 1)) count number of repetitions
        int numberOfRepetitions = 0;
        for (int index=0; index<MainTableSize; ++index) {
            if (table[index].raw()==rawRepeatedValue) { ++numberOfRepetitions; }
        }
        // 2)) print number of repetitions
        const int  expectedRepetitions = 1<<(8-data.length());
        const bool showRepetitions = numberOfRepetitions!=1 || expectedRepetitions!=1;
        if (showRepetitions) { std::printf("      %d/%d repetitions of %d,%d\n",
                                           numberOfRepetitions, expectedRepetitions, data.code(), data.length());
        }
        // 3)) print each repetition
        for (int index=0; index<MainTableSize; ++index) {
            if (table[index].raw()==rawRepeatedValue) {
                _printDebugInfo_TableEntry(table[index],index,0);
                isPrintedFlags[index] = true;
            }
        }
        if (showRepetitions) { std::printf("\n"); }
    }

    /**
     * Prints debugging information about a subtable of the internal huffman table
     */
    void ReversedHuffmanDecoder::_printDebugInfo_Subtable(Data* table, int index, bool* isPrintedFlags) {
        assert( table!=NULL && index<MainTableSize && isPrintedFlags!=NULL );
        assert( table[index]._isPointer() );
        const unsigned mask          = table[index]._mask();
        const int      firstSubindex = table[index]._index();
        const int      lastSubindex  = firstSubindex + mask;
        unsigned       bitstream     = index;
        
        _printDebugInfo_TableEntry( table[index], bitstream, mask );
        isPrintedFlags[index] = true;
        for (int subindex=firstSubindex; subindex<=lastSubindex; ++subindex) {
            _printDebugInfo_TableEntry( table[subindex], bitstream, mask );
            bitstream += 256;
        }
        std::printf("\n");
    }
    
    /**
     * Prints debugging information about an entry in the internal huffman table
     */
    void ReversedHuffmanDecoder::_printDebugInfo_TableEntry(Data data, unsigned bitStream, unsigned extraMask) {
        extraMask <<= 8;
        const bool     showExtraNumber = data.isValid();
        const bool     showSeparator   = data.isValid() ? data.length()>8 : true;
        const unsigned bracketPos      = data.isValid() && data.length()<=8 ? data.length() : 8;
        // 1)) print the bit sequence of the recognized bitstream
        for (int bitOrder=15; bitOrder>=8; --bitOrder) {
            std::printf("%c", !(extraMask&(1<<bitOrder)) ? ' ' :
                        !showExtraNumber           ? 'X' :
                        bitStream&(1<<bitOrder)    ? '1' : '0');
        }
        std::printf("%c", showSeparator ? '-' : ' ');
        for (int bitOrder=7; bitOrder>=0; --bitOrder) {
            if ((bitOrder+1)==bracketPos) { std::printf("[");  }
            std::printf("%c", bitStream&(1<<bitOrder) ? '1' : '0');
        }
        std::printf("]");
        
        // 2)) print the data attached to the recognized bitstream
        if      (data.raw()==0)  { std::printf(" <empty>\n");                              }
        else if (data.isValid()) { std::printf(" = %d,%d\n", data.code(), data.length());  }
        else                     { std::printf(" --> ptr to index [%d]\n", data._index()); }
    }


#   endif




}
