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
#ifndef PDZip_InvHuffmanDecoder_h
#define PDZip_InvHuffmanDecoder_h
#ifdef __cplusplus
#include <cassert>


typedef union InfHuff {
    struct value {
        unsigned length:15;
        unsigned isvalid:1;
        unsigned code:15;
    } value;
    struct subtable {
        unsigned mask:15;
        unsigned error:1;
        unsigned index:15;
    } subtable;
    unsigned raw;
} InfHuff;



namespace PDZip {


    /**
     * Represents a collection of code/length pairs
     *
     * Internally, the elements are kept sorted by length so it will be
     * more easy and efficient for ReverseHuffmanDecoder to process them.
     *
     * <...add more documentation...>
     *
     *
     * @date      Jul 26, 2018
     * @author    Martin Rizzo <martinrizzo@gmail.com>
     * @copyright Public Domain (CC0)
     */
    struct CodeLengths {
        
        static const int LastValidLength = 18;
        static const int LastValidCode   = 290;
        static const int TableSize      = (LastValidLength+1)+(LastValidCode+1);
        static const unsigned NextIndexMask = 0x03FF;
        static const unsigned CodeMask      = 0xFFFF;
        static const unsigned LengthMask    = 0x3F;
        
        //                         <  16bits  >:<6bits>:< 10bits  >
        // length,code,nextIndex = [   code   ]:[ len ]:[nextIndex]
        unsigned* _insertPtr[LastValidLength+1];
        unsigned  _table[TableSize];
        int       _nextIndex;
        int       _size;
        
        
        CodeLengths()
        : _nextIndex(0)
        { }

        inline bool isOpen() const {
            return (_nextIndex!=0);
        }

        inline void open() {
            assert( !isOpen() );
            for (int length=0; length<=LastValidLength; ++length) {
                *( _insertPtr[length] = &_table[length] ) = 0;
            }
            _size      = 0;
            _nextIndex = (LastValidLength+1);
            assert( isOpen() );
        }
        
        inline void close() {
            assert( isOpen() );
            int prevLength = 0;
            int length     = 0;
            do {
                unsigned nextIndex = 0;
                do {
                    nextIndex = (NextIndexMask & _table[++length]);
                } while ( length<LastValidLength && nextIndex==0 );
                *_insertPtr[prevLength] |= nextIndex;
                prevLength = length;
            } while ( length<LastValidLength );
            
            _nextIndex = 0;
            assert( !isOpen() );
        }
        
        inline void add(unsigned code, unsigned length) {
            assert( isOpen() );
            assert( 0<=code   &&   code<=LastValidCode   );
            assert( 0<=length && length<=LastValidLength );
            assert( _insertPtr[length]!=NULL );
            if (length>0) {
                *_insertPtr[length]  |= _nextIndex;
                *( _insertPtr[length] = &_table[ _nextIndex++ ] ) = code<<16 | length<<10;
            }
            ++_size;
        }
        
        inline void addFromInternalArray(unsigned char* internalArray, int numberOfElements) {
            assert( internalArray==getInternalArray(numberOfElements) );
            for (int i=0; i<numberOfElements; i++) {
                add(i, internalArray[i]);
            }
        }
        
        unsigned char* getInternalArray(int numberOfElements) {
            void* endOfArray = &_table[TableSize];
            return static_cast<unsigned char*>(endOfArray) - numberOfElements;
        }
        
        inline void addRange(int firstCode, int lastCode, unsigned length) {
            assert( isOpen() );
            assert( 0<=firstCode && firstCode<=(LastValidCode)   );
            assert( 0<=lastCode  &&  lastCode<=(LastValidCode+1) );
            assert( firstCode<=lastCode );
            for (int code=firstCode; code<lastCode; ++code) {
                add(code,length);
            }
        }
        
        inline const unsigned* getSourTable() const {
            return _table;
        }
        
        inline int size() const {
            return _size;
        }
        
        inline int getFirstIndex() const {
            assert( !isOpen() );
            return NextIndexMask & _table[0];
        }
        inline int getFirstIndexForLength(int length) const {
            assert( !isOpen() );
            assert( 0<=length && length<=LastValidLength );
            return NextIndexMask & _table[length];
        }
        inline int getNextIndex(int index) const {
            assert( !isOpen() );
            assert( 0<=index && index<TableSize );
            return NextIndexMask & _table[index];
        }
        
        inline unsigned getCodeAt(int index) const {
            assert( 0<=index && index<TableSize );
            return CodeMask & _table[index]>>16;
        }
        
        inline unsigned getLengthAt(int index) const {
            assert( 0<=index && index<TableSize );
            return LengthMask & _table[index]>>10;
        }
        
    public:
#   if !defined(NDEBUG)
        
        /**
         * Prints to console debug information about the current state of this object
         */
        void _printDebugInfo();
        
#   endif


    };
    

    /**
     * Huffman decoder for data with bits reversed
     *
     * <...add more documentation...>
     *
     *
     * @date      Jul 26, 2018
     * @author    Martin Rizzo <martinrizzo@gmail.com>
     * @copyright Public Domain (CC0)
     */
    class ReversedHuffmanDecoder {
        
        static const int MaxValidLength = 16;
        static const int MaxValidCode   = 32768;
        
        static const int MainTableSize = 256;    // < 8bits
        static const int FullTableSize = 2*1024; // < main-table + all sub-tables 

    public:
        /*
        class Data {
        public:
            Data() { }
            static Data fromCodeLength(unsigned code, unsigned length)  { return Data(length<<16 | 0x8000 | code); }
            static Data fromIndexMask(unsigned index, unsigned mask) { return Data(mask<<16 | index); }
            static Data fromRaw(unsigned raw) { return Data(raw); }
            bool     isValid()    const { return (_raw&0x8000); }
            unsigned code()       const { return (_raw&0x7FFF); }
            unsigned length()     const { return (_raw>>16);    }
            unsigned raw()        const { return _raw;          }
            bool     _isPointer() const { return !isValid(); }
            unsigned _index()     const { return (_raw&0x7FFF); }
            unsigned _mask()      const { return (_raw>>16);    }
        private:
            //                <   16bits   >:<   16bits   >
            // code,length  = 0000[ length ]:[1][  code   ]
            // reindex,mask = 0000[  mask  ]:[0][ reindex ]
            Data(unsigned raw) : _raw(raw) { }
            unsigned _raw;
        };
     */
        
    public:
        
        
        ReversedHuffmanDecoder() {
            _table[0].raw = 0;
        }
        
        bool isLoaded() const {
            return _table[0].raw != 0;
        }

        void load(const unsigned* table);

        /*
        void load(const CodeLengths& codeLengths);

         InfHuff decode8(unsigned bits8) const {
            return _table[bits8 & 0xFF];
        }
        
        InfHuff decode16(unsigned bits16, InfHuff data) const {
            assert( !data.subtable.error );
            return _table[ data.subtable.index + (bits16>>8 & data.subtable.mask) ];
        }

        InfHuff decode16(unsigned bits16) const {
            InfHuff data = _table[bits16 & 0xFF];
            return data.value.isvalid ? data : decode16(bits16,data);
        }
        */

    public:
#   if !defined(NDEBUG)
        
        /**
         * Prints to console the current state of the internal huffman table
         */
        void        printDebugInfo();
        static void _printDebugInfo_CodeRepetition(InfHuff* table, int firstIndex, bool* isPrintedFlags);
        static void _printDebugInfo_Subtable(InfHuff* table, int index, bool* isPrintedFlags);
        static void _printDebugInfo_TableEntry(InfHuff data, unsigned bitStream, unsigned extraMask);

#   endif

    public:
        InfHuff _table[FullTableSize];
    };
    
    

} /* namespace PDZip */
#endif /* ifdef __cplusplus */
#endif /* ifndef PDZip_InvHuffmanDecoder_h */
