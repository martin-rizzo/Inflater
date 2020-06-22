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

/* #define Data        PDZip::ReversedHuffmanDecoder::Data */
#define CodeLengths PDZip::CodeLengths


static inline unsigned reversedINC(unsigned huffman, unsigned length) {
    assert( length>0 );
    unsigned incBit = 1<<(length-1);
    do { huffman ^= incBit; } while ( (huffman&incBit)==0 && (incBit>>=1)!=0 );
    return huffman;
}

static inline void setWithRepetitions(InfHuff*  table,
                                      unsigned  availableEntries,
                                      unsigned  huffman,
                                      InfHuff   data,
                                      unsigned  length)
{
    assert( table!=NULL );
    assert( huffman < (1<<length) );
    assert( length>=1 );
    const unsigned unknownStep = (1<<length);
    for (unsigned unknownBits=0; unknownBits<availableEntries; unknownBits+=unknownStep) {
        table[ unknownBits | huffman ] = data;
    }
}

static inline int makeSubTable(InfHuff*           subtable,
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
        InfHuff data;
        data.value.code    = code;
        data.value.length  = length;
        data.value.isvalid = 1;
        setWithRepetitions(subtable, numberOfEntries, huffman>>DiscardedBits, data, length-DiscardedBits);

        huffman = reversedINC(huffman,length);
        index   = codeLengths.getNextIndex(index);
    }
    return numberOfEntries;
}

#undef Data
#undef CodeLengths


//==================================================================================================================
#pragma mark - CODE-LENGTHS RELATED FUNCTIONS

static const int      LastValidLength = 18;
static const int      LastValidCode   = 290;
static const int      TableSize       = (LastValidLength+1)+(LastValidCode+1);
static const unsigned NextIndexMask   = 0x03FF;
static const unsigned LengthMask      = 0x3F;
static const unsigned CodeMask        = 0xFFFF;


static int CL_getFirstIndex(const unsigned* sourtable) {
    assert( sourtable!=NULL );
    return NextIndexMask & sourtable[0];
}

static int CL_getNextIndex(const unsigned* sourtable, int index) {
    assert( sourtable!=NULL );
    assert( 0<=index && index<TableSize );
    return NextIndexMask & sourtable[index];
}


static unsigned CL_getCodeAt(const unsigned* sourtable, int index) {
    assert( sourtable!=NULL );
    assert( 0<=index && index<TableSize );
    return CodeMask & sourtable[index]>>16;
}

static unsigned CL_getLengthAt(const unsigned* sourtable, int index) {
    assert( sourtable!=NULL );
    assert( 0<=index && index<TableSize );
    return LengthMask & sourtable[index]>>10;
}


/* #define Data PDZip::ReversedHuffmanDecoder::Data */

static inline int CL_makeSubTable(InfHuff*           subtable,
                                  int                availableEntries,
                                  unsigned           firstHuffman,
                                  unsigned           maxLength,
                                  const unsigned*    sourtable,
                                  int                firstIndex,
                                  int                endIndex)
{
    static const unsigned DiscardedBits = 8; // < number of bits discarded by the subtable
    const int numberOfEntries = 1<<(maxLength-DiscardedBits);
    
    assert( numberOfEntries<=availableEntries  );
    
    unsigned huffman = firstHuffman;
    int      index   = firstIndex;
    while ( index!=endIndex ) {
        const unsigned code   = CL_getCodeAt(sourtable,index);
        const unsigned length = CL_getLengthAt(sourtable,index);
        InfHuff data;
        data.value.code    = code;
        data.value.length  = length;
        data.value.isvalid = 1;
        setWithRepetitions(subtable, numberOfEntries, huffman>>DiscardedBits, data, length-DiscardedBits);

        huffman = reversedINC(huffman,length);
        index   = CL_getNextIndex(sourtable,index);
    }
    return numberOfEntries;
}

#undef Data


namespace PDZip {

    
    //==================================================================================================================
#   pragma mark - LOADING THE HUFFMAN TABLE

void ReversedHuffmanDecoder::load(const unsigned int *sourtable) {
    assert( sourtable!=NULL  );
    
    // reset the main-table
    std::memset(_table, 0, MainTableSize*sizeof(InfHuff));
    
    unsigned huffman = 0;
    int      index   = CL_getFirstIndex(sourtable);
    unsigned length  = CL_getLengthAt(sourtable,index);

    // lengths from 1 to 8
    // unknown bits are filled with all possible values
    while ( index!=0 && length<=8 ) {
        InfHuff data;
        data.value.isvalid = 1;
        data.value.code    = CL_getCodeAt(sourtable,index);
        data.value.length  = length;
        setWithRepetitions(_table, 256, huffman, data, length);
        huffman = reversedINC(huffman,length);
        index   = CL_getNextIndex(sourtable,index);
        length  = CL_getLengthAt(sourtable,index);
    }
    // lengths from 9
    // subtables are created
    int insertIndex = MainTableSize;
    while ( index!=0 ) {
        const unsigned firstHuffman = huffman;
        const int      firstIndex   = index;
        do {
            length  = CL_getLengthAt(sourtable,index);
            huffman = reversedINC(huffman,length);
            index   = CL_getNextIndex(sourtable,index);
        } while ( index!=0 && (huffman&0xFF)==(firstHuffman&0xFF) );
        
        const unsigned maxLength = length;
        InfHuff data;
        data.subtable.error = 0;
        data.subtable.index = insertIndex;
        data.subtable.mask  = (1<<(maxLength-8))-1;
        _table[firstHuffman] = data;
        insertIndex += CL_makeSubTable(&_table[insertIndex],
                                       (FullTableSize-insertIndex),
                                       firstHuffman, maxLength,
                                       sourtable, firstIndex, index
                                       );
    }
}
    
    void ReversedHuffmanDecoder::load(const CodeLengths& codeLengths) {
        assert( codeLengths.isOpen()==false );
        
        // reset the main-table
        std::memset(_table, 0, MainTableSize*sizeof(InfHuff));
        
        unsigned huffman = 0;
        int      index   = codeLengths.getFirstIndex();
        unsigned length  = codeLengths.getLengthAt(index);

        // lengths from 1 to 8
        // unknown bits are filled with all possible values
        while ( index!=0 && length<=8 ) {
            InfHuff data;
            data.value.isvalid = 1;
            data.value.code    = codeLengths.getCodeAt(index);
            data.value.length  = length;
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
            InfHuff data;
            data.subtable.error = 0;
            data.subtable.index = insertIndex;
            data.subtable.mask  = (1<<(maxLength-8))-1;
            _table[firstHuffman] = data;
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
                if      ( _table[index].raw==0        ) { _printDebugInfo_TableEntry(_table[index], index, 0); isPrintedFlags[index]=true; }
                else if ( _table[index].value.isvalid ) { _printDebugInfo_CodeRepetition(_table, index, isPrintedFlags); }
                else                                    { _printDebugInfo_Subtable(_table, index, isPrintedFlags); }
            }
            index = index<255 ? reversedINC(index,8) : 256;
        }
        std::printf("\n");
    }

    /**
     * Prints debugging information about "code,length" repetitions in the internal huffman table
     */
    void ReversedHuffmanDecoder::_printDebugInfo_CodeRepetition(InfHuff* table, int firstIndex, bool* isPrintedFlags) {
        assert( table!=NULL && firstIndex<MainTableSize && isPrintedFlags!=NULL );
        assert( table[firstIndex].value.isvalid );
        const InfHuff  data             = table[firstIndex];
        const unsigned rawRepeatedValue = data.raw;
        
        // 1)) count number of repetitions
        int numberOfRepetitions = 0;
        for (int index=0; index<MainTableSize; ++index) {
            if (table[index].raw==rawRepeatedValue) { ++numberOfRepetitions; }
        }
        // 2)) print number of repetitions
        const int  expectedRepetitions = 1<<(8-data.value.length);
        const bool showRepetitions = numberOfRepetitions!=1 || expectedRepetitions!=1;
        if (showRepetitions) { std::printf("      %d/%d repetitions of %d,%d\n",
                                           numberOfRepetitions, expectedRepetitions, data.value.code, data.value.length);
        }
        // 3)) print each repetition
        for (int index=0; index<MainTableSize; ++index) {
            if (table[index].raw==rawRepeatedValue) {
                _printDebugInfo_TableEntry(table[index],index,0);
                isPrintedFlags[index] = true;
            }
        }
        if (showRepetitions) { std::printf("\n"); }
    }

    /**
     * Prints debugging information about a subtable of the internal huffman table
     */
    void ReversedHuffmanDecoder::_printDebugInfo_Subtable(InfHuff* table, int index, bool* isPrintedFlags) {
        assert( table!=NULL && index<MainTableSize && isPrintedFlags!=NULL );
        assert( !table[index].subtable.error  );
        const unsigned mask          = table[index].subtable.mask;
        const int      firstSubindex = table[index].subtable.index;
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
    void ReversedHuffmanDecoder::_printDebugInfo_TableEntry(InfHuff data, unsigned bitStream, unsigned extraMask) {
        extraMask <<= 8;
        const bool     showExtraNumber = data.value.isvalid;
        const bool     showSeparator   = data.value.isvalid ? data.value.length>8 : true;
        const unsigned bracketPos      = data.value.isvalid && data.value.length<=8 ? data.value.length : 8;
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
        if      (data.raw==0)        { std::printf(" <empty>\n");                                     }
        else if (data.value.isvalid) { std::printf(" = %d,%d\n", data.value.code, data.value.length); }
        else                         { std::printf(" --> ptr to index [%d]\n", data.subtable.index);  }
    }


#   endif




}
