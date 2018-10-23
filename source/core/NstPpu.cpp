////////////////////////////////////////////////////////////////////////////////////////
//
// Nestopia - NES/Famicom emulator written in C++
//
// Copyright (C) 2003-2008 Martin Freij
//
// This file is part of Nestopia.
//
// Nestopia is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Nestopia is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Nestopia; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cstring>

/*
decodePNG: The picoPNG function, decodes a PNG file buffer in memory, into a raw pixel buffer.
out_image: output parameter, this will contain the raw pixels after decoding.
  By default the output is 32-bit RGBA color.
  The std::vector is automatically resized to the correct size.
image_width: output_parameter, this will contain the width of the image in pixels.
image_height: output_parameter, this will contain the height of the image in pixels.
in_png: pointer to the buffer of the PNG file in memory. To get it from a file on
  disk, load it and store it in a memory buffer yourself first.
in_size: size of the input PNG file in bytes.
convert_to_rgba32: optional parameter, true by default.
  Set to true to get the output in RGBA 32-bit (8 bit per channel) color format
  no matter what color type the original PNG image had. This gives predictable,
  useable data from any random input PNG.
  Set to false to do no color conversion at all. The result then has the same data
  type as the PNG image, which can range from 1 bit to 64 bits per pixel.
  Information about the color type or palette colors are not provided. You need
  to know this information yourself to be able to use the data so this only
  works for trusted PNG files. Use LodePNG instead of picoPNG if you need this information.
return: 0 if success, not 0 if some error occured.
*/
namespace SkinHack
{
	int decodePNG(std::vector<unsigned char>& out_image, unsigned long& image_width, unsigned long& image_height, const unsigned char* in_png, size_t in_size, bool convert_to_rgba32 = true)
	{
	  // picoPNG version 20101224
	  // Copyright (c) 2005-2010 Lode Vandevenne
	  //
	  // This software is provided 'as-is', without any express or implied
	  // warranty. In no event will the authors be held liable for any damages
	  // arising from the use of this software.
	  //
	  // Permission is granted to anyone to use this software for any purpose,
	  // including commercial applications, and to alter it and redistribute it
	  // freely, subject to the following restrictions:
	  //
	  //     1. The origin of this software must not be misrepresented; you must not
	  //     claim that you wrote the original software. If you use this software
	  //     in a product, an acknowledgment in the product documentation would be
	  //     appreciated but is not required.
	  //     2. Altered source versions must be plainly marked as such, and must not be
	  //     misrepresented as being the original software.
	  //     3. This notice may not be removed or altered from any source distribution.
  
	  // picoPNG is a PNG decoder in one C++ function of around 500 lines. Use picoPNG for
	  // programs that need only 1 .cpp file. Since it's a single function, it's very limited,
	  // it can convert a PNG to raw pixel data either converted to 32-bit RGBA color or
	  // with no color conversion at all. For anything more complex, another tiny library
	  // is available: LodePNG (lodepng.c(pp)), which is a single source and header file.
	  // Apologies for the compact code style, it's to make this tiny.
  
	  static const unsigned long LENBASE[29] =  {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258};
	  static const unsigned long LENEXTRA[29] = {0,0,0,0,0,0,0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,  4,  5,  5,  5,  5,  0};
	  static const unsigned long DISTBASE[30] =  {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
	  static const unsigned long DISTEXTRA[30] = {0,0,0,0,1,1,2, 2, 3, 3, 4, 4, 5, 5,  6,  6,  7,  7,  8,  8,   9,   9,  10,  10,  11,  11,  12,   12,   13,   13};
	  static const unsigned long CLCL[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15}; //code length code lengths
	  struct Zlib //nested functions for zlib decompression
	  {
		static unsigned long readBitFromStream(size_t& bitp, const unsigned char* bits) { unsigned long result = (bits[bitp >> 3] >> (bitp & 0x7)) & 1; bitp++; return result;}
		static unsigned long readBitsFromStream(size_t& bitp, const unsigned char* bits, size_t nbits)
		{
		  unsigned long result = 0;
		  for(size_t i = 0; i < nbits; i++) result += (readBitFromStream(bitp, bits)) << i;
		  return result;
		}
		struct HuffmanTree
		{
		  int makeFromLengths(const std::vector<unsigned long>& bitlen, unsigned long maxbitlen)
		  { //make tree given the lengths
			unsigned long numcodes = (unsigned long)(bitlen.size()), treepos = 0, nodefilled = 0;
			std::vector<unsigned long> tree1d(numcodes), blcount(maxbitlen + 1, 0), nextcode(maxbitlen + 1, 0);
			for(unsigned long bits = 0; bits < numcodes; bits++) blcount[bitlen[bits]]++; //count number of instances of each code length
			for(unsigned long bits = 1; bits <= maxbitlen; bits++) nextcode[bits] = (nextcode[bits - 1] + blcount[bits - 1]) << 1;
			for(unsigned long n = 0; n < numcodes; n++) if(bitlen[n] != 0) tree1d[n] = nextcode[bitlen[n]]++; //generate all the codes
			tree2d.clear(); tree2d.resize(numcodes * 2, 32767); //32767 here means the tree2d isn't filled there yet
			for(unsigned long n = 0; n < numcodes; n++) //the codes
			for(unsigned long i = 0; i < bitlen[n]; i++) //the bits for this code
			{
			  unsigned long bit = (tree1d[n] >> (bitlen[n] - i - 1)) & 1;
			  if(treepos > numcodes - 2) return 55;
			  if(tree2d[2 * treepos + bit] == 32767) //not yet filled in
			  {
				if(i + 1 == bitlen[n]) { tree2d[2 * treepos + bit] = n; treepos = 0; } //last bit
				else { tree2d[2 * treepos + bit] = ++nodefilled + numcodes; treepos = nodefilled; } //addresses are encoded as values > numcodes
			  }
			  else treepos = tree2d[2 * treepos + bit] - numcodes; //subtract numcodes from address to get address value
			}
			return 0;
		  }
		  int decode(bool& decoded, unsigned long& result, size_t& treepos, unsigned long bit) const
		  { //Decodes a symbol from the tree
			unsigned long numcodes = (unsigned long)tree2d.size() / 2;
			if(treepos >= numcodes) return 11; //error: you appeared outside the codetree
			result = tree2d[2 * treepos + bit];
			decoded = (result < numcodes);
			treepos = decoded ? 0 : result - numcodes;
			return 0;
		  }
		  std::vector<unsigned long> tree2d; //2D representation of a huffman tree: The one dimension is "0" or "1", the other contains all nodes and leaves of the tree.
		};
		struct Inflator
		{
		  int error;
		  void inflate(std::vector<unsigned char>& out, const std::vector<unsigned char>& in, size_t inpos = 0)
		  {
			size_t bp = 0, pos = 0; //bit pointer and byte pointer
			error = 0;
			unsigned long BFINAL = 0;
			while(!BFINAL && !error)
			{
			  if(bp >> 3 >= in.size()) { error = 52; return; } //error, bit pointer will jump past memory
			  BFINAL = readBitFromStream(bp, &in[inpos]);
			  unsigned long BTYPE = readBitFromStream(bp, &in[inpos]); BTYPE += 2 * readBitFromStream(bp, &in[inpos]);
			  if(BTYPE == 3) { error = 20; return; } //error: invalid BTYPE
			  else if(BTYPE == 0) inflateNoCompression(out, &in[inpos], bp, pos, in.size());
			  else inflateHuffmanBlock(out, &in[inpos], bp, pos, in.size(), BTYPE);
			}
			if(!error) out.resize(pos); //Only now we know the true size of out, resize it to that
		  }
		  void generateFixedTrees(HuffmanTree& tree, HuffmanTree& treeD) //get the tree of a deflated block with fixed tree
		  {
			std::vector<unsigned long> bitlen(288, 8), bitlenD(32, 5);;
			for(size_t i = 144; i <= 255; i++) bitlen[i] = 9;
			for(size_t i = 256; i <= 279; i++) bitlen[i] = 7;
			tree.makeFromLengths(bitlen, 15);
			treeD.makeFromLengths(bitlenD, 15);
		  }
		  HuffmanTree codetree, codetreeD, codelengthcodetree; //the code tree for Huffman codes, dist codes, and code length codes
		  unsigned long huffmanDecodeSymbol(const unsigned char* in, size_t& bp, const HuffmanTree& codetree, size_t inlength)
		  { //decode a single symbol from given list of bits with given code tree. return value is the symbol
			bool decoded; unsigned long ct;
			for(size_t treepos = 0;;)
			{
			  if((bp & 0x07) == 0 && (bp >> 3) > inlength) { error = 10; return 0; } //error: end reached without endcode
			  error = codetree.decode(decoded, ct, treepos, readBitFromStream(bp, in)); if(error) return 0; //stop, an error happened
			  if(decoded) return ct;
			}
		  }
		  void getTreeInflateDynamic(HuffmanTree& tree, HuffmanTree& treeD, const unsigned char* in, size_t& bp, size_t inlength)
		  { //get the tree of a deflated block with dynamic tree, the tree itself is also Huffman compressed with a known tree
			std::vector<unsigned long> bitlen(288, 0), bitlenD(32, 0);
			if(bp >> 3 >= inlength - 2) { error = 49; return; } //the bit pointer is or will go past the memory
			size_t HLIT =  readBitsFromStream(bp, in, 5) + 257; //number of literal/length codes + 257
			size_t HDIST = readBitsFromStream(bp, in, 5) + 1; //number of dist codes + 1
			size_t HCLEN = readBitsFromStream(bp, in, 4) + 4; //number of code length codes + 4
			std::vector<unsigned long> codelengthcode(19); //lengths of tree to decode the lengths of the dynamic tree
			for(size_t i = 0; i < 19; i++) codelengthcode[CLCL[i]] = (i < HCLEN) ? readBitsFromStream(bp, in, 3) : 0;
			error = codelengthcodetree.makeFromLengths(codelengthcode, 7); if(error) return;
			size_t i = 0, replength;
			while(i < HLIT + HDIST)
			{
			  unsigned long code = huffmanDecodeSymbol(in, bp, codelengthcodetree, inlength); if(error) return;
			  if(code <= 15)  { if(i < HLIT) bitlen[i++] = code; else bitlenD[i++ - HLIT] = code; } //a length code
			  else if(code == 16) //repeat previous
			  {
				if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
				replength = 3 + readBitsFromStream(bp, in, 2);
				unsigned long value; //set value to the previous code
				if((i - 1) < HLIT) value = bitlen[i - 1];
				else value = bitlenD[i - HLIT - 1];
				for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
				{
				  if(i >= HLIT + HDIST) { error = 13; return; } //error: i is larger than the amount of codes
				  if(i < HLIT) bitlen[i++] = value; else bitlenD[i++ - HLIT] = value;
				}
			  }
			  else if(code == 17) //repeat "0" 3-10 times
			  {
				if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
				replength = 3 + readBitsFromStream(bp, in, 3);
				for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
				{
				  if(i >= HLIT + HDIST) { error = 14; return; } //error: i is larger than the amount of codes
				  if(i < HLIT) bitlen[i++] = 0; else bitlenD[i++ - HLIT] = 0;
				}
			  }
			  else if(code == 18) //repeat "0" 11-138 times
			  {
				if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
				replength = 11 + readBitsFromStream(bp, in, 7);
				for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
				{
				  if(i >= HLIT + HDIST) { error = 15; return; } //error: i is larger than the amount of codes
				  if(i < HLIT) bitlen[i++] = 0; else bitlenD[i++ - HLIT] = 0;
				}
			  }
			  else { error = 16; return; } //error: somehow an unexisting code appeared. This can never happen.
			}
			if(bitlen[256] == 0) { error = 64; return; } //the length of the end code 256 must be larger than 0
			error = tree.makeFromLengths(bitlen, 15); if(error) return; //now we've finally got HLIT and HDIST, so generate the code trees, and the function is done
			error = treeD.makeFromLengths(bitlenD, 15); if(error) return;
		  }
		  void inflateHuffmanBlock(std::vector<unsigned char>& out, const unsigned char* in, size_t& bp, size_t& pos, size_t inlength, unsigned long btype) 
		  {
			if(btype == 1) { generateFixedTrees(codetree, codetreeD); }
			else if(btype == 2) { getTreeInflateDynamic(codetree, codetreeD, in, bp, inlength); if(error) return; }
			for(;;)
			{
			  unsigned long code = huffmanDecodeSymbol(in, bp, codetree, inlength); if(error) return;
			  if(code == 256) return; //end code
			  else if(code <= 255) //literal symbol
			  {
				if(pos >= out.size()) out.resize((pos + 1) * 2); //reserve more room
				out[pos++] = (unsigned char)(code);
			  }
			  else if(code >= 257 && code <= 285) //length code
			  {
				size_t length = LENBASE[code - 257], numextrabits = LENEXTRA[code - 257];
				if((bp >> 3) >= inlength) { error = 51; return; } //error, bit pointer will jump past memory
				length += readBitsFromStream(bp, in, numextrabits);
				unsigned long codeD = huffmanDecodeSymbol(in, bp, codetreeD, inlength); if(error) return;
				if(codeD > 29) { error = 18; return; } //error: invalid dist code (30-31 are never used)
				unsigned long dist = DISTBASE[codeD], numextrabitsD = DISTEXTRA[codeD];
				if((bp >> 3) >= inlength) { error = 51; return; } //error, bit pointer will jump past memory
				dist += readBitsFromStream(bp, in, numextrabitsD);
				size_t start = pos, back = start - dist; //backwards
				if(pos + length >= out.size()) out.resize((pos + length) * 2); //reserve more room
				for(size_t i = 0; i < length; i++) { out[pos++] = out[back++]; if(back >= start) back = start - dist; }
			  }
			}
		  }
		  void inflateNoCompression(std::vector<unsigned char>& out, const unsigned char* in, size_t& bp, size_t& pos, size_t inlength)
		  {
			while((bp & 0x7) != 0) bp++; //go to first boundary of byte
			size_t p = bp / 8;
			if(p >= inlength - 4) { error = 52; return; } //error, bit pointer will jump past memory
			unsigned long LEN = in[p] + 256 * in[p + 1], NLEN = in[p + 2] + 256 * in[p + 3]; p += 4;
			if(LEN + NLEN != 65535) { error = 21; return; } //error: NLEN is not one's complement of LEN
			if(pos + LEN >= out.size()) out.resize(pos + LEN);
			if(p + LEN > inlength) { error = 23; return; } //error: reading outside of in buffer
			for(unsigned long n = 0; n < LEN; n++) out[pos++] = in[p++]; //read LEN bytes of literal data
			bp = p * 8;
		  }
		};
		int decompress(std::vector<unsigned char>& out, const std::vector<unsigned char>& in) //returns error value
		{
		  Inflator inflator;
		  if(in.size() < 2) { return 53; } //error, size of zlib data too small
		  if((in[0] * 256 + in[1]) % 31 != 0) { return 24; } //error: 256 * in[0] + in[1] must be a multiple of 31, the FCHECK value is supposed to be made that way
		  unsigned long CM = in[0] & 15, CINFO = (in[0] >> 4) & 15, FDICT = (in[1] >> 5) & 1;
		  if(CM != 8 || CINFO > 7) { return 25; } //error: only compression method 8: inflate with sliding window of 32k is supported by the PNG spec
		  if(FDICT != 0) { return 26; } //error: the specification of PNG says about the zlib stream: "The additional flags shall not specify a preset dictionary."
		  inflator.inflate(out, in, 2);
		  return inflator.error; //note: adler32 checksum was skipped and ignored
		}
	  };
	  struct PNG //nested functions for PNG decoding
	  {
		struct Info
		{
		  unsigned long width, height, colorType, bitDepth, compressionMethod, filterMethod, interlaceMethod, key_r, key_g, key_b;
		  bool key_defined; //is a transparent color key given?
		  std::vector<unsigned char> palette;
		} info;
		int error;
		void decode(std::vector<unsigned char>& out, const unsigned char* in, size_t size, bool convert_to_rgba32)
		{
		  error = 0;
		  if(size == 0 || in == 0) { error = 48; return; } //the given data is empty
		  readPngHeader(&in[0], size); if(error) return;
		  size_t pos = 33; //first byte of the first chunk after the header
		  std::vector<unsigned char> idat; //the data from idat chunks
		  bool IEND = false, known_type = true;
		  info.key_defined = false;
		  while(!IEND) //loop through the chunks, ignoring unknown chunks and stopping at IEND chunk. IDAT data is put at the start of the in buffer
		  {
			if(pos + 8 >= size) { error = 30; return; } //error: size of the in buffer too small to contain next chunk
			size_t chunkLength = read32bitInt(&in[pos]); pos += 4;
			if(chunkLength > 2147483647) { error = 63; return; }
			if(pos + chunkLength >= size) { error = 35; return; } //error: size of the in buffer too small to contain next chunk
			if(in[pos + 0] == 'I' && in[pos + 1] == 'D' && in[pos + 2] == 'A' && in[pos + 3] == 'T') //IDAT chunk, containing compressed image data
			{
			  idat.insert(idat.end(), &in[pos + 4], &in[pos + 4 + chunkLength]);
			  pos += (4 + chunkLength);
			}
			else if(in[pos + 0] == 'I' && in[pos + 1] == 'E' && in[pos + 2] == 'N' && in[pos + 3] == 'D')  { pos += 4; IEND = true; }
			else if(in[pos + 0] == 'P' && in[pos + 1] == 'L' && in[pos + 2] == 'T' && in[pos + 3] == 'E') //palette chunk (PLTE)
			{
			  pos += 4; //go after the 4 letters
			  info.palette.resize(4 * (chunkLength / 3));
			  if(info.palette.size() > (4 * 256)) { error = 38; return; } //error: palette too big
			  for(size_t i = 0; i < info.palette.size(); i += 4)
			  {
				for(size_t j = 0; j < 3; j++) info.palette[i + j] = in[pos++]; //RGB
				info.palette[i + 3] = 255; //alpha
			  }
			}
			else if(in[pos + 0] == 't' && in[pos + 1] == 'R' && in[pos + 2] == 'N' && in[pos + 3] == 'S') //palette transparency chunk (tRNS)
			{
			  pos += 4; //go after the 4 letters
			  if(info.colorType == 3)
			  {
				if(4 * chunkLength > info.palette.size()) { error = 39; return; } //error: more alpha values given than there are palette entries
				for(size_t i = 0; i < chunkLength; i++) info.palette[4 * i + 3] = in[pos++];
			  }
			  else if(info.colorType == 0)
			  {
				if(chunkLength != 2) { error = 40; return; } //error: this chunk must be 2 bytes for greyscale image
				info.key_defined = 1; info.key_r = info.key_g = info.key_b = 256 * in[pos] + in[pos + 1]; pos += 2;
			  }
			  else if(info.colorType == 2)
			  {
				if(chunkLength != 6) { error = 41; return; } //error: this chunk must be 6 bytes for RGB image
				info.key_defined = 1;
				info.key_r = 256 * in[pos] + in[pos + 1]; pos += 2;
				info.key_g = 256 * in[pos] + in[pos + 1]; pos += 2;
				info.key_b = 256 * in[pos] + in[pos + 1]; pos += 2;
			  }
			  else { error = 42; return; } //error: tRNS chunk not allowed for other color models
			}
			else //it's not an implemented chunk type, so ignore it: skip over the data
			{
			  if(!(in[pos + 0] & 32)) { error = 69; return; } //error: unknown critical chunk (5th bit of first byte of chunk type is 0)
			  pos += (chunkLength + 4); //skip 4 letters and uninterpreted data of unimplemented chunk
			  known_type = false;
			}
			pos += 4; //step over CRC (which is ignored)
		  }
		  unsigned long bpp = getBpp(info);
		  std::vector<unsigned char> scanlines(((info.width * (info.height * bpp + 7)) / 8) + info.height); //now the out buffer will be filled
		  Zlib zlib; //decompress with the Zlib decompressor
		  error = zlib.decompress(scanlines, idat); if(error) return; //stop if the zlib decompressor returned an error
		  size_t bytewidth = (bpp + 7) / 8, outlength = (info.height * info.width * bpp + 7) / 8;
		  out.resize(outlength); //time to fill the out buffer
		  unsigned char* out_ = outlength ? &out[0] : 0; //use a regular pointer to the std::vector for faster code if compiled without optimization
		  if(info.interlaceMethod == 0) //no interlace, just filter
		  {
			size_t linestart = 0, linelength = (info.width * bpp + 7) / 8; //length in bytes of a scanline, excluding the filtertype byte
			if(bpp >= 8) //byte per byte
			for(unsigned long y = 0; y < info.height; y++)
			{
			  unsigned long filterType = scanlines[linestart];
			  const unsigned char* prevline = (y == 0) ? 0 : &out_[(y - 1) * info.width * bytewidth];
			  unFilterScanline(&out_[linestart - y], &scanlines[linestart + 1], prevline, bytewidth, filterType,  linelength); if(error) return;
			  linestart += (1 + linelength); //go to start of next scanline
			}
			else //less than 8 bits per pixel, so fill it up bit per bit
			{
			  std::vector<unsigned char> templine((info.width * bpp + 7) >> 3); //only used if bpp < 8
			  for(size_t y = 0, obp = 0; y < info.height; y++)
			  {
				unsigned long filterType = scanlines[linestart];
				const unsigned char* prevline = (y == 0) ? 0 : &out_[(y - 1) * info.width * bytewidth];
				unFilterScanline(&templine[0], &scanlines[linestart + 1], prevline, bytewidth, filterType, linelength); if(error) return;
				for(size_t bp = 0; bp < info.width * bpp;) setBitOfReversedStream(obp, out_, readBitFromReversedStream(bp, &templine[0]));
				linestart += (1 + linelength); //go to start of next scanline
			  }
			}
		  }
		  else //interlaceMethod is 1 (Adam7)
		  {
			size_t passw[7] = { (info.width + 7) / 8, (info.width + 3) / 8, (info.width + 3) / 4, (info.width + 1) / 4, (info.width + 1) / 2, (info.width + 0) / 2, (info.width + 0) / 1 };
			size_t passh[7] = { (info.height + 7) / 8, (info.height + 7) / 8, (info.height + 3) / 8, (info.height + 3) / 4, (info.height + 1) / 4, (info.height + 1) / 2, (info.height + 0) / 2 };
			size_t passstart[7] = {0};
			size_t pattern[28] = {0,4,0,2,0,1,0,0,0,4,0,2,0,1,8,8,4,4,2,2,1,8,8,8,4,4,2,2}; //values for the adam7 passes
			for(int i = 0; i < 6; i++) passstart[i + 1] = passstart[i] + passh[i] * ((passw[i] ? 1 : 0) + (passw[i] * bpp + 7) / 8);
			std::vector<unsigned char> scanlineo((info.width * bpp + 7) / 8), scanlinen((info.width * bpp + 7) / 8); //"old" and "new" scanline
			for(int i = 0; i < 7; i++)
			  adam7Pass(&out_[0], &scanlinen[0], &scanlineo[0], &scanlines[passstart[i]], info.width, pattern[i], pattern[i + 7], pattern[i + 14], pattern[i + 21], passw[i], passh[i], bpp);
		  }
		  if(convert_to_rgba32 && (info.colorType != 6 || info.bitDepth != 8)) //conversion needed
		  {
			std::vector<unsigned char> data = out;
			error = convert(out, &data[0], info, info.width, info.height);
		  }
		}
		void readPngHeader(const unsigned char* in, size_t inlength) //read the information from the header and store it in the Info
		{
		  if(inlength < 29) { error = 27; return; } //error: the data length is smaller than the length of the header
		  if(in[0] != 137 || in[1] != 80 || in[2] != 78 || in[3] != 71 || in[4] != 13 || in[5] != 10 || in[6] != 26 || in[7] != 10) { error = 28; return; } //no PNG signature
		  if(in[12] != 'I' || in[13] != 'H' || in[14] != 'D' || in[15] != 'R') { error = 29; return; } //error: it doesn't start with a IHDR chunk!
		  info.width = read32bitInt(&in[16]); info.height = read32bitInt(&in[20]);
		  info.bitDepth = in[24]; info.colorType = in[25];
		  info.compressionMethod = in[26]; if(in[26] != 0) { error = 32; return; } //error: only compression method 0 is allowed in the specification
		  info.filterMethod = in[27]; if(in[27] != 0) { error = 33; return; } //error: only filter method 0 is allowed in the specification
		  info.interlaceMethod = in[28]; if(in[28] > 1) { error = 34; return; } //error: only interlace methods 0 and 1 exist in the specification
		  error = checkColorValidity(info.colorType, info.bitDepth);
		}
		void unFilterScanline(unsigned char* recon, const unsigned char* scanline, const unsigned char* precon, size_t bytewidth, unsigned long filterType, size_t length)
		{
		  switch(filterType)
		  {
			case 0: for(size_t i = 0; i < length; i++) recon[i] = scanline[i]; break;
			case 1:
			  for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i];
			  for(size_t i = bytewidth; i <    length; i++) recon[i] = (scanline[i] + recon[i - bytewidth]) & 0xFF;
			  break;
			case 2:
			  if(precon) for(size_t i = 0; i < length; i++) recon[i] = (scanline[i] + precon[i]) & 0xFF;
			  else       for(size_t i = 0; i < length; i++) recon[i] = scanline[i];
			  break;
			case 3:
			  if(precon)
			  {
				for(size_t i =         0; i < bytewidth; i++) recon[i] = (scanline[i] + precon[i] / 2) & 0xFF;
				for(size_t i = bytewidth; i <    length; i++) recon[i] = (scanline[i] + ((recon[i - bytewidth] + precon[i]) / 2)) & 0xFF;
			  }
			  else
			  {
				for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i];
				for(size_t i = bytewidth; i <    length; i++) recon[i] = (scanline[i] + recon[i - bytewidth] / 2) & 0xFF;
			  }
			  break;
			case 4:
			  if(precon)
			  {
				for(size_t i =         0; i < bytewidth; i++) recon[i] = (scanline[i] + paethPredictor(0, precon[i], 0)) & 0xFF;
				for(size_t i = bytewidth; i <    length; i++) recon[i] = (scanline[i] + paethPredictor(recon[i - bytewidth], precon[i], precon[i - bytewidth])) & 0xFF;
			  }
			  else
			  {
				for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i];
				for(size_t i = bytewidth; i <    length; i++) recon[i] = (scanline[i] + paethPredictor(recon[i - bytewidth], 0, 0)) & 0xFF;
			  }
			  break;
			default: error = 36; return; //error: unexisting filter type given
		  }
		}
		void adam7Pass(unsigned char* out, unsigned char* linen, unsigned char* lineo, const unsigned char* in, unsigned long w, size_t passleft, size_t passtop, size_t spacex, size_t spacey, size_t passw, size_t passh, unsigned long bpp)
		{ //filter and reposition the pixels into the output when the image is Adam7 interlaced. This function can only do it after the full image is already decoded. The out buffer must have the correct allocated memory size already.
		  if(passw == 0) return;
		  size_t bytewidth = (bpp + 7) / 8, linelength = 1 + ((bpp * passw + 7) / 8);
		  for(unsigned long y = 0; y < passh; y++)
		  {
			unsigned char filterType = in[y * linelength], *prevline = (y == 0) ? 0 : lineo;
			unFilterScanline(linen, &in[y * linelength + 1], prevline, bytewidth, filterType, (w * bpp + 7) / 8); if(error) return;
			if(bpp >= 8) for(size_t i = 0; i < passw; i++) for(size_t b = 0; b < bytewidth; b++) //b = current byte of this pixel
			  out[bytewidth * w * (passtop + spacey * y) + bytewidth * (passleft + spacex * i) + b] = linen[bytewidth * i + b];
			else for(size_t i = 0; i < passw; i++)
			{
			  size_t obp = bpp * w * (passtop + spacey * y) + bpp * (passleft + spacex * i), bp = i * bpp;
			  for(size_t b = 0; b < bpp; b++) setBitOfReversedStream(obp, out, readBitFromReversedStream(bp, &linen[0]));
			}
			unsigned char* temp = linen; linen = lineo; lineo = temp; //swap the two buffer pointers "line old" and "line new"
		  }
		}
		static unsigned long readBitFromReversedStream(size_t& bitp, const unsigned char* bits) { unsigned long result = (bits[bitp >> 3] >> (7 - (bitp & 0x7))) & 1; bitp++; return result;}
		static unsigned long readBitsFromReversedStream(size_t& bitp, const unsigned char* bits, unsigned long nbits)
		{
		  unsigned long result = 0;
		  for(size_t i = nbits - 1; i < nbits; i--) result += ((readBitFromReversedStream(bitp, bits)) << i);
		  return result;
		}
		void setBitOfReversedStream(size_t& bitp, unsigned char* bits, unsigned long bit) { bits[bitp >> 3] |=  (bit << (7 - (bitp & 0x7))); bitp++; }
		unsigned long read32bitInt(const unsigned char* buffer) { return (buffer[0] << 24) | (buffer[1] << 16) | (buffer[2] << 8) | buffer[3]; }
		int checkColorValidity(unsigned long colorType, unsigned long bd) //return type is a LodePNG error code
		{
		  if((colorType == 2 || colorType == 4 || colorType == 6)) { if(!(bd == 8 || bd == 16)) return 37; else return 0; }
		  else if(colorType == 0) { if(!(bd == 1 || bd == 2 || bd == 4 || bd == 8 || bd == 16)) return 37; else return 0; }
		  else if(colorType == 3) { if(!(bd == 1 || bd == 2 || bd == 4 || bd == 8            )) return 37; else return 0; }
		  else return 31; //unexisting color type
		}
		unsigned long getBpp(const Info& info)
		{
		  if(info.colorType == 2) return (3 * info.bitDepth);
		  else if(info.colorType >= 4) return (info.colorType - 2) * info.bitDepth;
		  else return info.bitDepth;
		}
		int convert(std::vector<unsigned char>& out, const unsigned char* in, Info& infoIn, unsigned long w, unsigned long h)
		{ //converts from any color type to 32-bit. return value = LodePNG error code
		  size_t numpixels = w * h, bp = 0;
		  out.resize(numpixels * 4);
		  unsigned char* out_ = out.empty() ? 0 : &out[0]; //faster if compiled without optimization
		  if(infoIn.bitDepth == 8 && infoIn.colorType == 0) //greyscale
		  for(size_t i = 0; i < numpixels; i++)
		  {
			out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[i];
			out_[4 * i + 3] = (infoIn.key_defined && in[i] == infoIn.key_r) ? 0 : 255;
		  }
		  else if(infoIn.bitDepth == 8 && infoIn.colorType == 2) //RGB color
		  for(size_t i = 0; i < numpixels; i++)
		  {
			for(size_t c = 0; c < 3; c++) out_[4 * i + c] = in[3 * i + c];
			out_[4 * i + 3] = (infoIn.key_defined == 1 && in[3 * i + 0] == infoIn.key_r && in[3 * i + 1] == infoIn.key_g && in[3 * i + 2] == infoIn.key_b) ? 0 : 255;
		  }
		  else if(infoIn.bitDepth == 8 && infoIn.colorType == 3) //indexed color (palette)
		  for(size_t i = 0; i < numpixels; i++)
		  {
			if(4U * in[i] >= infoIn.palette.size()) return 46;
			for(size_t c = 0; c < 4; c++) out_[4 * i + c] = infoIn.palette[4 * in[i] + c]; //get rgb colors from the palette
		  }
		  else if(infoIn.bitDepth == 8 && infoIn.colorType == 4) //greyscale with alpha
		  for(size_t i = 0; i < numpixels; i++)
		  {
			out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[2 * i + 0];
			out_[4 * i + 3] = in[2 * i + 1];
		  }
		  else if(infoIn.bitDepth == 8 && infoIn.colorType == 6) for(size_t i = 0; i < numpixels; i++) for(size_t c = 0; c < 4; c++) out_[4 * i + c] = in[4 * i + c]; //RGB with alpha
		  else if(infoIn.bitDepth == 16 && infoIn.colorType == 0) //greyscale
		  for(size_t i = 0; i < numpixels; i++)
		  {
			out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[2 * i];
			out_[4 * i + 3] = (infoIn.key_defined && 256U * in[i] + in[i + 1] == infoIn.key_r) ? 0 : 255;
		  }
		  else if(infoIn.bitDepth == 16 && infoIn.colorType == 2) //RGB color
		  for(size_t i = 0; i < numpixels; i++)
		  {
			for(size_t c = 0; c < 3; c++) out_[4 * i + c] = in[6 * i + 2 * c];
			out_[4 * i + 3] = (infoIn.key_defined && 256U*in[6*i+0]+in[6*i+1] == infoIn.key_r && 256U*in[6*i+2]+in[6*i+3] == infoIn.key_g && 256U*in[6*i+4]+in[6*i+5] == infoIn.key_b) ? 0 : 255;
		  }
		  else if(infoIn.bitDepth == 16 && infoIn.colorType == 4) //greyscale with alpha
		  for(size_t i = 0; i < numpixels; i++)
		  {
			out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[4 * i]; //most significant byte
			out_[4 * i + 3] = in[4 * i + 2];
		  }
		  else if(infoIn.bitDepth == 16 && infoIn.colorType == 6) for(size_t i = 0; i < numpixels; i++) for(size_t c = 0; c < 4; c++) out_[4 * i + c] = in[8 * i + 2 * c]; //RGB with alpha
		  else if(infoIn.bitDepth < 8 && infoIn.colorType == 0) //greyscale
		  for(size_t i = 0; i < numpixels; i++)
		  {
			unsigned long value = (readBitsFromReversedStream(bp, in, infoIn.bitDepth) * 255) / ((1 << infoIn.bitDepth) - 1); //scale value from 0 to 255
			out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = (unsigned char)(value);
			out_[4 * i + 3] = (infoIn.key_defined && value && ((1U << infoIn.bitDepth) - 1U) == infoIn.key_r && ((1U << infoIn.bitDepth) - 1U)) ? 0 : 255;
		  }
		  else if(infoIn.bitDepth < 8 && infoIn.colorType == 3) //palette
		  for(size_t i = 0; i < numpixels; i++)
		  {
			unsigned long value = readBitsFromReversedStream(bp, in, infoIn.bitDepth);
			if(4 * value >= infoIn.palette.size()) return 47;
			for(size_t c = 0; c < 4; c++) out_[4 * i + c] = infoIn.palette[4 * value + c]; //get rgb colors from the palette
		  }
		  return 0;
		}
		unsigned char paethPredictor(short a, short b, short c) //Paeth predicter, used by PNG filter type 4
		{
		  short p = a + b - c, pa = p > a ? (p - a) : (a - p), pb = p > b ? (p - b) : (b - p), pc = p > c ? (p - c) : (c - p);
		  return (unsigned char)((pa <= pb && pa <= pc) ? a : pb <= pc ? b : c);
		}
	  };
	  PNG decoder; decoder.decode(out_image, in_png, in_size, convert_to_rgba32);
	  image_width = decoder.info.width; image_height = decoder.info.height;
	  return decoder.error;
	}

	void loadFile(std::vector<unsigned char>& buffer, const std::string& filename) //designed for loading files from hard disk in an std::vector
	{
	  std::ifstream file(filename.c_str(), std::ios::in|std::ios::binary|std::ios::ate);

	  //get filesize
	  std::streamsize size = 0;
	  if(file.seekg(0, std::ios::end).good()) size = file.tellg();
	  if(file.seekg(0, std::ios::beg).good()) size -= file.tellg();

	  //read contents of the file into the vector
	  if(size > 0)
	  {
		buffer.resize((size_t)size);
		file.read((char*)(&buffer[0]), size);
	  }
	  else buffer.clear();
	}

	struct ReplacementPattern {
		unsigned short& operator[](unsigned int idx) { return data[idx]; }
		unsigned short data[8];
		unsigned char  isTrigger;
	};

	std::map<unsigned int,ReplacementPattern> patternLut[128*16*2];

	bool enable = true;

	unsigned char getNesColor(unsigned int rgba)
	{
		switch(rgba & 0x00FCFCFC) { // clean up transparency info and 2 lowest "service" bits of each color
			case 0x747474:	return 0x0;
			case 0x8c1824:	return 0x1;
			case 0xa80000:	return 0x2;
			case 0x9c0044:	return 0x3;
			case 0x74008c:	return 0x4;
			case 0x1000a8:	return 0x5;
			case 0x0000a4:	return 0x6;
			case 0x00087c:	return 0x7;
			case 0x002c40:	return 0x8;
			case 0x004400:	return 0x9;
			case 0x005000:	return 0xa;
			case 0x143c00:	return 0xb;
			case 0x5c3c18:	return 0xc;
			case 0xbcbcbc:	return 0x10;
			case 0xec7000:	return 0x11;
			case 0xec3820:	return 0x12;
			case 0xf00080:	return 0x13;
			case 0xbc00bc:	return 0x14;
			case 0x5800e4:	return 0x15;
			case 0x0028d8:	return 0x16;
			case 0x0c4cc8:	return 0x17;
			case 0x007088:	return 0x18;
			case 0x009400:	return 0x19;
			case 0x00a800:	return 0x1a;
			case 0x389000:	return 0x1b;
			case 0x888000:	return 0x1c;
			case 0xfcbc3c:	return 0x21;
			case 0xfc945c:	return 0x22;
			case 0xfc88cc:	return 0x23;
			case 0xfc78f4:	return 0x24;
			case 0xb474fc:	return 0x25;
			case 0x6074fc:	return 0x26;
			case 0x3898fc:	return 0x27;
			case 0x3cbcf0:	return 0x28;
			case 0x10d080:	return 0x29;
			case 0x48dc4c:	return 0x2a;
			case 0x98f858:	return 0x2b;
			case 0xd8e800:	return 0x2c;
			case 0xfce4a8:	return 0x31;
			case 0xfcd4c4:	return 0x32;
			case 0xfcc8d4:	return 0x33;
			case 0xfcc4fc:	return 0x34;
			case 0xd8c4fc:	return 0x35;
			case 0xb0bcfc:	return 0x36;
			case 0xa8d8fc:	return 0x37;
			case 0xa0e4fc:	return 0x38;
			case 0xa0fce0:	return 0x39;
			case 0xbcf0a8:	return 0x3a;
			case 0xccfcb0:	return 0x3b;
			case 0xf0fc9c:	return 0x3c;
			default:		return 0;
		}
	}

	void parseSkin()
	{
		for( int k=0; k<4096; k++ ) {
			patternLut[k].clear();
		}
		char filename[] = "skin.png";
		std::ofstream log;
		log.open ("skin.log");
		log << std::hex;
		log << "Trying to parse " << filename << std::endl;

		//load and decode
		std::vector<unsigned char> buffer, image;
		loadFile(buffer, filename);

		if( buffer.empty() ) {
			log << "File not found or empty; Quit parsing" << std::endl;
			log.close();
			return;
		}

		unsigned long w, h;
		int error = decodePNG(image, w, h, buffer.empty() ? 0 : &buffer[0], (unsigned long)buffer.size());
  
		//if there's an error, log it
		if(error != 0) {
			log << "Error loading PNG: " << error << std::endl;
			buffer.clear();
			log.close();
			return;
		}
  
		//the pixels are now in the vector "image", use it as texture, draw it, ...
		if(w != 256 || (h & 0xFF)) {
			log << "Error: PNG should be 256 by 256*x pixels" << std::endl;
			buffer.clear();
			image.clear();
			log.close();
			return;
		}

		unsigned char trigger = 0;
		bool await_trigger = false;
		// number of 256 pixel tall skins
		int num_skins = h >> 8;
		for( int s=0; s<num_skins; s++ ) {
			log << "Parsing skin N: " << s << std::endl;
			unsigned int *rgba = (unsigned int*)&image[s*256*256*4];
			// converting original pattern images into original patterns 
			// each original pattern image contains 8 32bit rgba values
			// each original pattern is a dword, where lsw is 8 pixels, 2 bits (4 colors) each
			// two rightmost bits of a msw contain an original pattern attribute

			// also converting replacement pattern images into replacement patterns
			// each replacement pattern image contains 8 32bit rgba values
			// each replacement pattern contains 8 16bit values, where msb==1 for non-transparent colors, then go 5 bits red, 5 bits green and 5 bits blue

			// if first original pattern is not defined (its first pixel is transparent), and the trigger was not set in a previous skin, 
			// then we await a trigger pattern in a current skin.
			// first defined original pattern will be considered the trigger pattern for the current skin (and next skins if they also lack pattern 0)
			// other patterns of the current skin will be replaced only if the trigger pattern was met before when rendering a screen
			await_trigger = !(rgba[0] & 0xFF000000);
			if (!await_trigger) {
				trigger = 0; // reset the trigger, if the current skin has pattern 0 filled
			}

			for( int p=0; p<4096*8; p+=8 ) {
				if (!(rgba[p] & 0xFF000000)) { // original pattern not defined (its first pixel is transparent)
					continue;
				}
				bool replacement_transparent = true;
				// pattern attribute is stored in two rightmost bits of a last green byte of an original pattern image
				unsigned int original_pattern = (rgba[p + 7] >> 8) & 3;
				ReplacementPattern replacement_pattern;
				for( int x=0; x<8; x++ ) {
					original_pattern <<= 2;
					// pattern color is stored in two rightmost bits of a red byte of an original pattern image pixel
					original_pattern |= rgba[p + x] & 3;
					// replacement pixels are stored 128 rows lower than original pixels
					unsigned int replacement_pixel = rgba[128*256 + p + x];
					// each replacement pattern contains 8 16bit values, where msb==1, then go 5 bits red, 5 bits green and 5 bits blue
					if (!(replacement_pixel & 0xFF000000)) { // replacement pattern image pixel is considered transparent if its alpha equals 0
						replacement_pattern[x] = 0;
					} else {
						replacement_transparent = false;
						replacement_pattern[x] = 0x8000 
																+ ((replacement_pixel & 0x0000F8) << 7)  // b
																+ ((replacement_pixel & 0x00F800) >> 6)  // g
																+ ((replacement_pixel & 0xF80000) >> 19); // r
					}
				}
				// log << "Patt.N: " << (p>>3) << std::endl;
				// if replacement pattern is not transparent, store it in a replacement pattern lut with the original pattern as a key
				// first, calculate a correct offset 
				// patterns in the same column of rows 0-7,8-15... should follow each other, so that p==0 has the offset 0, 8->8, 16->16... 256->1, 264->9...
				// int offset = p%256 + (p%2048)/256 + (p/2048)*256;
				int offset = (p & 0xFF) + ((p & 0x7FF) >> 8) + ((p & ~0x7FF) >> 3);
				// leftmost 16 patterns in every row (of 32) go to bank 0 with offset 0, rightmost 16 patterns go to bank 1 with offset 0x800
				// offset = (offset%256 > 127)*2048 + (offset/256)*128 + offset%128;
				offset = ((offset & 0x80) << 4) + ((offset & ~0xFF) >> 1) + (offset & 0x7F);
				// then insert original and replacement patterns as a key-value pair into the correct offset of the lut
				// log << "Offset: " << offset << std::endl;
				if (!replacement_transparent) {
					// Nestopia for its convinience puts odd pixels into msb and even pixels into lsb
					unsigned int nestopia_pattern = 
						(original_pattern & 0x30000) |
						((original_pattern & 0xc000) >> 8) |
						((original_pattern & 0x3000) << 2) |
						((original_pattern & 0x0c00) >> 6) |
						((original_pattern & 0x0300) << 4) |
						((original_pattern & 0x00c0) >> 4) |
						((original_pattern & 0x0030) << 6) |
						((original_pattern & 0x000c) >> 2) |
						((original_pattern & 0x0003) << 8);

					// Add a NES color index of the 6th pixel to a lookup pattern; This is to fix a problem with incorrect replacement of solid fills
					// By now it's clear that I should've just constructed a lookup pattern from 8 6bit NES color indices
					unsigned char last_pixel_color = getNesColor(rgba[p + 6]);
					nestopia_pattern |= last_pixel_color << 18;

					nestopia_pattern |= trigger << 24;
					if (!trigger && await_trigger) {
						trigger = (offset >> 3) & 0xFF;
						// log << "Trigger: " << (int)trigger << std::endl;
						replacement_pattern.isTrigger = (rgba[128*256 + p] & 0xFF000000) < 0xFF000000 ? 2 : 1; // isTrigger == 2 - a trigger that won't work if hClock is more than 247; Castlevania stage 10 hack
						// log << "Trigger type: " << (int)replacement_pattern.isTrigger << std::endl;
					} else {
						replacement_pattern.isTrigger = 0;
					}

					patternLut[offset].insert(std::make_pair(nestopia_pattern,replacement_pattern));
					// log << "Original: " << original_pattern << std::endl;
					// log << "Nestopia: " << nestopia_pattern << std::endl;
					// log << "Last clr: " << (int)last_pixel_color << std::endl;
					// log << "Replacement: ";
					for( int i=0; i<8; i++ ){
						// log << replacement_pattern[i] << " ";
					}
					// log << std::endl;
				}
			}
			log << "Done." << std::endl;
		}
		log.close();
		buffer.clear();
		image.clear();
	}

	unsigned char triggerPattern;
	void clearTriggerPattern()
	{
		triggerPattern = 0;
	}

	bool tryFetchBg(unsigned int address, unsigned char pattern0, unsigned char pattern1, unsigned char attribute, unsigned short* bgPixels, unsigned short* palette, unsigned int hclock)
	{
		if (address & 8) {
			return false; // do not try to fetch from an upper bit plane
		}
		
		unsigned char last_pixel_color = (unsigned char)palette[(pattern0 & 3) | ((attribute & 3) << 2)];
		if( (last_pixel_color & 0xf) > 0xc || last_pixel_color == 0x20 || last_pixel_color == 0x30 ) {
			last_pixel_color = 0;
		}
		unsigned int pattern_lut_offset = ((address & ~7) >> 1) + (address & 7);
		unsigned int original_pattern = (triggerPattern << 24) | (last_pixel_color << 18) | ((attribute & 3) << 16) | (pattern1 << 8) | pattern0;
		std::map<unsigned int,ReplacementPattern>::iterator it = patternLut[pattern_lut_offset].find(original_pattern);
		if( it == patternLut[pattern_lut_offset].end() ) {
			if( !triggerPattern ) {
				return false;
			} else {
				original_pattern &= 0xFFFFFF; // if a triggered pattern is not found, try to find a regular replacement pattern
				it = patternLut[pattern_lut_offset].find(original_pattern);
			}
		} 
		if( it == patternLut[pattern_lut_offset].end() ) {
			return false;
		} else {
			std::memcpy(bgPixels, &(it->second[0]), sizeof(it->second.data));
			// found a trigger
			if( it->second.isTrigger && (hclock < 248 || it->second.isTrigger == 1)) {
				triggerPattern = (pattern_lut_offset >> 3) & 0xFF; // set a trigger until the end of a frame
			}
			return true;
		}
	}

	bool tryFetchSprite(unsigned int address, unsigned int pattern, unsigned char attribute, unsigned char flip, unsigned short* spritePixels, unsigned short* palette)
	{
		unsigned char last_pixel_color = (pattern & 3) ? (unsigned char)palette[16 + ((pattern & 3) | ((attribute & 3) << 2))] : 0;
		if( (last_pixel_color & 0xf) > 0xc || last_pixel_color == 0x20 || last_pixel_color == 0x30 ) {
			last_pixel_color = 0;
		}
		unsigned int pattern_lut_offset = ((address & ~15) >> 1) + (address & 7);
		unsigned int original_pattern = (triggerPattern << 24) | (last_pixel_color << 18) | ((attribute & 3) << 16) | pattern;
		std::map<unsigned int,ReplacementPattern>::iterator it = patternLut[pattern_lut_offset].find(original_pattern);
		if( it == patternLut[pattern_lut_offset].end() ) {
			if( !triggerPattern ) {
				return false;
			} else {
				original_pattern &= 0xFFFFFF; // if a triggered pattern is not found, try to find a regular replacement pattern
				it = patternLut[pattern_lut_offset].find(original_pattern);
			}
		} 
		if( it == patternLut[pattern_lut_offset].end() ) {
			return false;
		} else {
			if( flip ) {
				unsigned short *ptr = &(it->second[0]);
				for( int i=0; i<8; i++ ) {
					spritePixels[i] = ptr[7-i];
				}
			}
			else {
				std::memcpy(spritePixels, &(it->second[0]), sizeof(it->second.data));
			}
			return true;
		}
	}
}

#include "NstCpu.hpp"
#include "NstPpu.hpp"
#include "NstState.hpp"

namespace Nes
{
	namespace Core
	{
		#ifdef NST_MSVC_OPTIMIZE
		#pragma optimize("s", on)
		#endif

		const byte Ppu::yuvMaps[4][0x40] =
		{
			{
				0x35, 0x23, 0x16, 0x22, 0x1C, 0x09, 0x2D, 0x15,
				0x20, 0x00, 0x27, 0x05, 0x04, 0x28, 0x08, 0x20,
				0x21, 0x27, 0x07, 0x29, 0x3C, 0x32, 0x36, 0x12,
				0x28, 0x2B, 0x0D, 0x08, 0x10, 0x3D, 0x24, 0x01,
				0x01, 0x31, 0x33, 0x2A, 0x2C, 0x0C, 0x1B, 0x14,
				0x0D, 0x07, 0x34, 0x06, 0x13, 0x02, 0x26, 0x0D,
				0x0D, 0x19, 0x10, 0x0A, 0x39, 0x03, 0x37, 0x17,
				0x09, 0x11, 0x1A, 0x1D, 0x38, 0x25, 0x18, 0x3A
			},
			{
				0x0D, 0x27, 0x18, 0x39, 0x3A, 0x25, 0x1C, 0x31,
				0x16, 0x13, 0x38, 0x34, 0x20, 0x23, 0x3C, 0x1A,
				0x09, 0x21, 0x06, 0x10, 0x1B, 0x29, 0x08, 0x22,
				0x2D, 0x24, 0x01, 0x2B, 0x32, 0x08, 0x0D, 0x03,
				0x04, 0x36, 0x26, 0x33, 0x11, 0x07, 0x10, 0x02,
				0x14, 0x28, 0x00, 0x09, 0x12, 0x0D, 0x28, 0x20,
				0x27, 0x1D, 0x2A, 0x17, 0x0C, 0x01, 0x15, 0x19,
				0x0D, 0x2C, 0x07, 0x37, 0x35, 0x05, 0x0A, 0x3D
			},
			{
				0x14, 0x25, 0x3A, 0x10, 0x1A, 0x20, 0x31, 0x09,
				0x01, 0x0D, 0x36, 0x08, 0x15, 0x10, 0x27, 0x3C,
				0x22, 0x1C, 0x05, 0x12, 0x19, 0x18, 0x17, 0x1B,
				0x00, 0x03, 0x0D, 0x02, 0x16, 0x06, 0x34, 0x35,
				0x23, 0x09, 0x01, 0x37, 0x1D, 0x27, 0x26, 0x20,
				0x29, 0x04, 0x21, 0x24, 0x11, 0x3D, 0x0D, 0x07,
				0x2C, 0x08, 0x39, 0x33, 0x07, 0x2A, 0x28, 0x2D,
				0x0A, 0x0D, 0x32, 0x38, 0x13, 0x2B, 0x28, 0x0C
			},
			{
				0x18, 0x03, 0x1C, 0x28, 0x0D, 0x35, 0x01, 0x17,
				0x10, 0x07, 0x2A, 0x01, 0x36, 0x37, 0x1A, 0x39,
				0x25, 0x08, 0x12, 0x34, 0x0D, 0x2D, 0x06, 0x26,
				0x27, 0x1B, 0x22, 0x19, 0x04, 0x0D, 0x3A, 0x21,
				0x05, 0x0A, 0x07, 0x02, 0x13, 0x14, 0x00, 0x15,
				0x0C, 0x10, 0x11, 0x09, 0x1D, 0x38, 0x3D, 0x24,
				0x33, 0x20, 0x08, 0x16, 0x28, 0x2B, 0x20, 0x3C,
				0x0D, 0x27, 0x23, 0x31, 0x29, 0x32, 0x2C, 0x09
			}
		};

		Ppu::Tiles::Tiles()
		: padding0(0), padding1(0) {}

		Ppu::Oam::Oam()
		: limit(buffer + STD_LINE_SPRITES*4), spriteLimit(true) {}

		Ppu::Output::Output(Video::Screen::Pixel* p)
		: pixels(p) {}

		Ppu::TileLut::TileLut()
		{
			for (uint i=0; i < 0x400; ++i)
			{
				block[i][0] = (i & 0xC0) ? (i >> 6 & 0xC) | (i >> 6 & 0x3) : 0;
				block[i][1] = (i & 0x30) ? (i >> 6 & 0xC) | (i >> 4 & 0x3) : 0;
				block[i][2] = (i & 0x0C) ? (i >> 6 & 0xC) | (i >> 2 & 0x3) : 0;
				block[i][3] = (i & 0x03) ? (i >> 6 & 0xC) | (i >> 0 & 0x3) : 0;
			}
		}

		Ppu::Ppu(Cpu& c)
		:
		cpu    (c),
		output (screen.pixels),
		model  (PPU_RP2C02),
		rgbMap (NULL),
		yuvMap (NULL)
		{
			cycles.one = PPU_RP2C02_CC;
			PowerOff();
			SkinHack::parseSkin();
		}

		void Ppu::PowerOff()
		{
			Reset( true, false, false );
		}

		void Ppu::Reset(bool hard,bool acknowledged)
		{
			Reset( hard, acknowledged, true );
		}

		void Ppu::Reset(const bool hard,const bool acknowledged,const bool map)
		{
			if (map)
			{
				for (uint i=0x2000; i < 0x4000; i += 0x8)
				{
					cpu.Map( i+0 ).Set( this, i != 0x3000 ? &Ppu::Peek_2xxx : &Ppu::Peek_3000, &Ppu::Poke_2000 );
					cpu.Map( i+1 ).Set( this,               &Ppu::Peek_2xxx,                   &Ppu::Poke_2001 );
					cpu.Map( i+2 ).Set( this,               &Ppu::Peek_2002,                   &Ppu::Poke_2xxx );
					cpu.Map( i+3 ).Set( this,               &Ppu::Peek_2xxx,                   &Ppu::Poke_2003 );
					cpu.Map( i+4 ).Set( this,               &Ppu::Peek_2004,                   &Ppu::Poke_2004 );
					cpu.Map( i+5 ).Set( this,               &Ppu::Peek_2xxx,                   &Ppu::Poke_2005 );
					cpu.Map( i+6 ).Set( this,               &Ppu::Peek_2xxx,                   &Ppu::Poke_2006 );
					cpu.Map( i+7 ).Set( this,               &Ppu::Peek_2007,                   &Ppu::Poke_2007 );
				}

				switch (model)
				{
					case PPU_RC2C05_01:
					case PPU_RC2C05_02:
					case PPU_RC2C05_03:
					case PPU_RC2C05_04:

						if (model == PPU_RC2C05_02)
						{
							for (uint i=0x2002; i < 0x4000; i += 0x8)
								cpu.Map( i ).Set( &Ppu::Peek_2002_RC2C05_02 );
						}
						else if (model == PPU_RC2C05_03)
						{
							for (uint i=0x2002; i < 0x4000; i += 0x8)
								cpu.Map( i ).Set( &Ppu::Peek_2002_RC2C05_03 );
						}
						else
						{
							for (uint i=0x2002; i < 0x4000; i += 0x8)
								cpu.Map( i ).Set( &Ppu::Peek_2002_RC2C05_01_04 );
						}

					case PPU_RC2C05_05:

						for (uint i=0x2000; i < 0x4000; i += 0x8)
						{
							cpu.Map( i+0 ).Set( &Ppu::Poke_2001 );
							cpu.Map( i+1 ).Set( &Ppu::Poke_2000 );
						}
						break;
				}

				cpu.Map( 0x4014U ).Set( this, &Ppu::Peek_4014, &Ppu::Poke_4014 );
			}

			if (hard)
			{
				static const byte powerUpPalette[] =
				{
					0x3F,0x01,0x00,0x01,0x00,0x02,0x02,0x0D,
					0x08,0x10,0x08,0x24,0x00,0x00,0x04,0x2C,
					0x09,0x01,0x34,0x03,0x00,0x04,0x00,0x14,
					0x08,0x3A,0x00,0x02,0x00,0x20,0x2C,0x08
				};

				std::memcpy( palette.ram, powerUpPalette, Palette::SIZE );
				std::memset( oam.ram, Oam::GARBAGE, Oam::SIZE );
				std::memset( nameTable.ram, NameTable::GARBAGE, NameTable::SIZE );

				io.latch = 0;
				io.buffer = Io::BUFFER_GARBAGE;

				regs.status = 0;
				regs.ctrl[0] = 0;
				regs.ctrl[1] = 0;
				regs.frame = 0;
				regs.oam = 0;

				scroll.latch = 0;
				scroll.xFine = 0;
				scroll.toggle = 0;
				scroll.address = 0;

				output.burstPhase = 0;

				cycles.reset = 0;
				cycles.hClock = HCLOCK_BOOT;
			}
			else if (acknowledged)
			{
				io.buffer = 0;

				regs.status = 0;
				regs.ctrl[0] = 0;
				regs.ctrl[1] = 0;

				scroll.latch = 0;
				scroll.xFine = 0;
				scroll.toggle = 0;

				cycles.reset = Cpu::CYCLE_MAX;
				cycles.hClock = HCLOCK_BOOT;

				std::memset( oam.ram, Oam::GARBAGE, Oam::SIZE );
			}
			else
			{
				cycles.hClock = HCLOCK_DUMMY;
				cycles.reset = 0;
			}

			if (chr.Source().Empty())
			{
				chr.Source().Set( Ram::RAM, true, false, NameTable::SIZE, nameTable.ram );
				chr.SwapBanks<SIZE_2K,0x0000>(0,0,0,0);
			}

			if (nmt.Source().Empty())
			{
				nmt.Source().Set( Ram::RAM, true, true, NameTable::SIZE, nameTable.ram );
				nmt.SwapBanks<SIZE_2K,0x0000>(0,0);
			}

			chr.ResetAccessor();
			nmt.ResetAccessors();

			cycles.vClock = 0;
			cycles.count = Cpu::CYCLE_MAX;

			scanline = SCANLINE_VBLANK;

			io.address = 0;
			io.pattern = 0;
			io.line.Unset();

			tiles.pattern[0] = 0;
			tiles.pattern[1] = 0;
			tiles.attribute = 0;
			tiles.index = 8;
			tiles.mask = 0;

			oam.index = 0;
			oam.address = 0;
			oam.latch = 0;
			oam.spriteZeroInLine = false;
			oam.phase = &Ppu::EvaluateSpritesPhase0;
			oam.buffered = oam.buffer;
			oam.visible = oam.output;
			oam.mask = 0;

			output.target = NULL;

			hActiveHook.Unset();
			hBlankHook.Unset();

			UpdateStates();

			screen.Clear();
		}

		uint Ppu::SetAddressLineHook(const Core::Io::Line& line)
		{
			io.line = line;
			return io.address;
		}

		void Ppu::SetHActiveHook(const Hook& hook)
		{
			hActiveHook = hook;
		}

		void Ppu::SetHBlankHook(const Hook& hook)
		{
			hBlankHook = hook;
		}

		void Ppu::UpdateStates()
		{
			oam.height = (regs.ctrl[0] >> 2 & 8) + 8;

			/*original code begin*/
			//tiles.show[0] = (regs.ctrl[1] & Regs::CTRL1_BG_ENABLED) ? 0xFF : 0x00;
			//tiles.show[1] = (regs.ctrl[1] & Regs::CTRL1_BG_ENABLED_NO_CLIP) == Regs::CTRL1_BG_ENABLED_NO_CLIP ? 0xFF : 0x00;

			//oam.show[0] = (regs.ctrl[1] & Regs::CTRL1_SP_ENABLED) ? 0xFF : 0x00;
			//oam.show[1] = (regs.ctrl[1] & Regs::CTRL1_SP_ENABLED_NO_CLIP) == Regs::CTRL1_SP_ENABLED_NO_CLIP ? 0xFF : 0x00;
			/*original code end*/
			/*skin hack begin*/
			tiles.show[0] = (regs.ctrl[1] & Regs::CTRL1_BG_ENABLED) ? 0xFFFF : 0x00;
			tiles.show[1] = (regs.ctrl[1] & Regs::CTRL1_BG_ENABLED_NO_CLIP) == Regs::CTRL1_BG_ENABLED_NO_CLIP ? 0xFFFF : 0x00;

			oam.show[0] = (regs.ctrl[1] & Regs::CTRL1_SP_ENABLED) ? 0xFFFF : 0x00;
			oam.show[1] = (regs.ctrl[1] & Regs::CTRL1_SP_ENABLED_NO_CLIP) == Regs::CTRL1_SP_ENABLED_NO_CLIP ? 0xFFFF : 0x00;
			/*skin hack end*/

			UpdatePalette();
		}

		void Ppu::UpdatePalette()
		{
			for (uint i=0, c=Coloring(), e=Emphasis(); i < Palette::SIZE; ++i)
				output.palette[i] = (rgbMap ? rgbMap[palette.ram[i] & uint(Palette::COLOR)] : palette.ram[i]) & c | e;
		}

		void Ppu::SaveState(State::Saver& state,const dword baseChunk) const
		{
			state.Begin( baseChunk );

			{
				const byte data[11] =
				{
					regs.ctrl[0],
					regs.ctrl[1],
					regs.status,
					scroll.address & 0xFF,
					scroll.address >> 8,
					scroll.latch & 0xFF,
					scroll.latch >> 8,
					scroll.xFine | scroll.toggle << 3,
					regs.oam,
					io.buffer,
					io.latch
				};

				state.Begin( AsciiId<'R','E','G'>::V ).Write( data ).End();
			}

			state.Begin( AsciiId<'P','A','L'>::V ).Compress( palette.ram   ).End();
			state.Begin( AsciiId<'O','A','M'>::V ).Compress( oam.ram       ).End();
			state.Begin( AsciiId<'N','M','T'>::V ).Compress( nameTable.ram ).End();

			if (model == PPU_RP2C02)
				state.Begin( AsciiId<'F','R','M'>::V ).Write8( (regs.frame & Regs::FRAME_ODD) == 0 ).End();

			if (cycles.hClock == HCLOCK_BOOT)
				state.Begin( AsciiId<'P','O','W'>::V ).Write8( 0x0 ).End();

			state.End();
		}

		void Ppu::LoadState(State::Loader& state)
		{
			cycles.hClock = HCLOCK_DUMMY;
			regs.frame = 0;
			output.burstPhase = 0;

			while (const dword chunk = state.Begin())
			{
				switch (chunk)
				{
					case AsciiId<'R','E','G'>::V:
					{
						State::Loader::Data<11> data( state );

						regs.ctrl[0]   = data[0];
						regs.ctrl[1]   = data[1];
						regs.status    = data[2] & Regs::STATUS_BITS;
						scroll.address = data[3] | (data[4] << 8 & 0x7F00);
						scroll.latch   = data[5] | (data[6] << 8 & 0x7F00);
						scroll.xFine   = data[7] & 0x7;
						scroll.toggle  = data[7] >> 3 & 0x1;
						regs.oam       = data[8];
						io.buffer      = data[9];
						io.latch       = data[10];

						break;
					}

					case AsciiId<'P','A','L'>::V:

						state.Uncompress( palette.ram );
						break;

					case AsciiId<'O','A','M'>::V:

						state.Uncompress( oam.ram );
						break;

					case AsciiId<'N','M','T'>::V:

						state.Uncompress( nameTable.ram );
						break;

					case AsciiId<'F','R','M'>::V:

						if (model == PPU_RP2C02)
							regs.frame = (state.Read8() & 0x1) ? 0 : Regs::FRAME_ODD;

						break;

					case AsciiId<'P','O','W'>::V:

						cycles.hClock = HCLOCK_BOOT;
						break;
				}

				state.End();
			}

			// SkinHack::parseSkin();

			UpdateStates();
		}

		void Ppu::EnableCpuSynchronization()
		{
			cpu.AddHook( Hook(this,&Ppu::Hook_Sync) );
		}

		void Ppu::ChrMem::ResetAccessor()
		{
			accessor.Set( this, &ChrMem::Access_Pattern );
		}

		void Ppu::NmtMem::ResetAccessors()
		{
			accessors[0].Set( this, &NmtMem::Access_Name_2000 );
			accessors[1].Set( this, &NmtMem::Access_Name_2400 );
			accessors[2].Set( this, &NmtMem::Access_Name_2800 );
			accessors[3].Set( this, &NmtMem::Access_Name_2C00 );
		}

		void Ppu::SetModel(const PpuModel m,const bool yuvConversion)
		{
			if (model != m)
			{
				model = m;
				regs.frame = 0;
				output.burstPhase = 0;

				switch (model)
				{
					case PPU_RP2C07: cycles.one = PPU_RP2C07_CC; break;
					case PPU_DENDY:  cycles.one = PPU_DENDY_CC;  break;
					default:         cycles.one = PPU_RP2C02_CC; break;
				}
			}

			const byte* const map =
			(
				model == PPU_RP2C04_0001 ? yuvMaps[0] :
				model == PPU_RP2C04_0002 ? yuvMaps[1] :
				model == PPU_RP2C04_0003 ? yuvMaps[2] :
				model == PPU_RP2C04_0004 ? yuvMaps[3] :
                                           NULL
			);

			const byte* const tmp[2] =
			{
				yuvConversion ? NULL : map,
				yuvConversion ? map : NULL
			};

			if (yuvMap != tmp[0] || rgbMap != tmp[1])
			{
				yuvMap = tmp[0];
				rgbMap = tmp[1];

				UpdatePalette();
			}
		}

		#ifdef NST_MSVC_OPTIMIZE
		#pragma optimize("", on)
		#endif

		NST_FORCE_INLINE Cycle Ppu::GetCycles() const
		{
			return (cycles.vClock + cycles.hClock) * cycles.one;
		}

		NST_FORCE_INLINE Cycle Ppu::GetLocalCycles(Cycle clock) const
		{
			NST_COMPILE_ASSERT( PPU_DENDY_CC == PPU_RP2C02_CC || PPU_DENDY_CC == PPU_RP2C07_CC );
			return cycles.one == PPU_RP2C02_CC ? clock / PPU_RP2C02_CC : (clock+PPU_RP2C07_CC-1) / PPU_RP2C07_CC;
		}

		void Ppu::BeginFrame(bool frameLock)
		{
			NST_ASSERT
			(
				(scanline == SCANLINE_VBLANK) &&
				(cycles.hClock == HCLOCK_BOOT || cycles.hClock == HCLOCK_DUMMY) &&
				(cpu.GetModel() == CPU_RP2A07) == (model == PPU_RP2C07) &&
				(cpu.GetModel() == CPU_DENDY)  == (model == PPU_DENDY)
			);

			oam.limit = oam.buffer + ((oam.spriteLimit || frameLock) ? Oam::STD_LINE_SPRITES*4 : Oam::MAX_LINE_SPRITES*4);
			output.target = output.pixels;

			Cycle frame;

			switch (model)
			{
				case PPU_RP2C02:

					regs.frame ^= Regs::FRAME_ODD;

				default:

					if (cycles.hClock == HCLOCK_DUMMY)
					{
						cycles.vClock = PPU_RP2C02_HVINT / PPU_RP2C02_CC - HCLOCK_DUMMY;
						cycles.count = PPU_RP2C02_HVINT;
						frame = PPU_RP2C02_HVSYNC_0;
					}
					else
					{
						cycles.vClock = PPU_RP2C02_HVSYNCBOOT / PPU_RP2C02_CC - HCLOCK_BOOT;
						cycles.count = PPU_RP2C02_HVSYNCBOOT;
						frame = PPU_RP2C02_HVSYNCBOOT;
					}
					break;

				case PPU_RP2C07:

					if (cycles.hClock == HCLOCK_DUMMY)
					{
						cycles.vClock = PPU_RP2C07_HVINT / PPU_RP2C07_CC - HCLOCK_DUMMY;
						cycles.count = PPU_RP2C07_HVINT;
						frame = PPU_RP2C07_HVSYNC;
					}
					else
					{
						cycles.vClock = PPU_RP2C07_HVSYNCBOOT / PPU_RP2C07_CC - HCLOCK_BOOT;
						cycles.count = PPU_RP2C07_HVSYNCBOOT;
						frame = PPU_RP2C07_HVSYNCBOOT;
					}
					break;

				case PPU_DENDY:

					if (cycles.hClock == HCLOCK_DUMMY)
					{
						cycles.vClock = PPU_DENDY_HVINT / PPU_DENDY_CC - HCLOCK_DUMMY;
						cycles.count = PPU_DENDY_HVINT;
						frame = PPU_DENDY_HVSYNC;
					}
					else
					{
						cycles.vClock = PPU_DENDY_HVSYNCBOOT / PPU_DENDY_CC - HCLOCK_BOOT;
						cycles.count = PPU_DENDY_HVSYNCBOOT;
						frame = PPU_DENDY_HVSYNCBOOT;
					}
					break;
			}

			cpu.SetFrameCycles( frame );
			SkinHack::clearTriggerPattern();
		}

		NES_HOOK(Ppu,Sync)
		{
			const Cycle elapsed = cpu.GetCycles();

			if (cycles.count < elapsed)
			{
				cycles.count = GetLocalCycles( elapsed ) - cycles.vClock;
				Run();
			}
		}

		void Ppu::EndFrame()
		{
			if (cycles.count != Cpu::CYCLE_MAX)
			{
				cycles.count = Cpu::CYCLE_MAX;
				Run();
			}
		}

		void Ppu::Update(Cycle dataSetup,const uint readAddress)
		{
			dataSetup += cpu.Update( readAddress );

			if (cycles.count < dataSetup)
			{
				cycles.count = GetLocalCycles( dataSetup ) - cycles.vClock;
				Run();
			}
		}

		void Ppu::SetMirroring(const byte (&banks)[4])
		{
			Update( cycles.one );

			NST_ASSERT( banks[0] < 4 && banks[1] < 4 && banks[2] < 4 && banks[3] < 4 );
			nmt.SwapBanks<SIZE_1K,0x0000>( banks[0], banks[1], banks[2], banks[3] );
		}

		void Ppu::SetMirroring(NmtMirroring mirroring)
		{
			Update( cycles.one );

			nmt.SwapBanks<SIZE_1K,0x0000>
			(
				uint(mirroring) >> 0 & 0x1U,
				uint(mirroring) >> 1 & 0x1U,
				uint(mirroring) >> 2 & 0x1U,
				uint(mirroring) >> 3 & 0x1U
			);
		}

		NES_ACCESSOR(Ppu::ChrMem,Pattern)
		{
			return Peek( address );
		}

		NES_ACCESSOR(Ppu::NmtMem,Name_2000)
		{
			return (*this)[0][address];
		}

		NES_ACCESSOR(Ppu::NmtMem,Name_2400)
		{
			return (*this)[1][address];
		}

		NES_ACCESSOR(Ppu::NmtMem,Name_2800)
		{
			return (*this)[2][address];
		}

		NES_ACCESSOR(Ppu::NmtMem,Name_2C00)
		{
			return (*this)[3][address];
		}

		NST_FORCE_INLINE uint Ppu::Chr::FetchPattern(uint address) const
		{
			return accessor.Fetch( address & 0x1FFF );
		}

		NST_FORCE_INLINE uint Ppu::Nmt::FetchName(uint address) const
		{
			return accessors[address >> 10 & 0x3].Fetch( address & 0x3FF );
		}

		NST_FORCE_INLINE uint Ppu::Nmt::FetchAttribute(uint address) const
		{
			return accessors[address >> 10 & 0x3].Fetch( 0x3C0 | (address & 0x03F) );
		}

		NST_FORCE_INLINE void Ppu::UpdateAddressLine(uint address)
		{
			NST_ASSERT( address <= 0x3FFF );
			io.address = address;

			if (io.line)
				io.line.Toggle( io.address, GetCycles() );
		}

		NST_FORCE_INLINE void Ppu::UpdateVramAddress()
		{
			if ((scanline != SCANLINE_VBLANK ) && (regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED))
			{
				scroll.ClockX();
				scroll.ClockY();
			}
			else
			{
				scroll.address = (scroll.address + ((regs.ctrl[0] & Regs::CTRL0_INC32) ? 32 : 1)) & 0x7FFF;
			}
		}

		NST_FORCE_INLINE void Ppu::OpenName()
		{
			UpdateAddressLine( 0x2000 | (scroll.address & 0x0FFF) );
		}

		NST_FORCE_INLINE void Ppu::FetchName()
		{
			io.pattern = nmt.FetchName( io.address ) << 4 | scroll.address >> 12 | (regs.ctrl[0] << 8 & 0x1000);
		}

		NST_FORCE_INLINE void Ppu::OpenAttribute()
		{
			UpdateAddressLine( 0x23C0 | (scroll.address & 0x0C00) | (scroll.address >> 4 & 0x0038) | (scroll.address >> 2 & 0x0007) );
		}

		NST_FORCE_INLINE void Ppu::FetchAttribute()
		{
			tiles.attribute = nmt.FetchAttribute( io.address ) >> ((scroll.address & 0x2) | (scroll.address >> 4 & 0x4));
		}

		NST_FORCE_INLINE void Ppu::OpenPattern(uint address)
		{
			UpdateAddressLine( address );
		}

		NST_FORCE_INLINE uint Ppu::FetchSpPattern() const
		{
			return chr.FetchPattern( io.address );
		}

		NST_FORCE_INLINE void Ppu::FetchBgPattern0()
		{
			const uint pattern = chr.FetchPattern( io.address );

			tiles.pattern[1] = pattern >> 0 & 0x55;
			tiles.pattern[0] = pattern >> 1 & 0x55;
		}

		NST_FORCE_INLINE void Ppu::FetchBgPattern1()
		{
			const uint pattern = chr.FetchPattern( io.address );

			tiles.pattern[0] |= pattern << 0 & 0xAA;
			tiles.pattern[1] |= pattern << 1 & 0xAA;
		}

		uint Ppu::GetPixelCycles() const
		{
			return (scanline+1)-1U < 240 ? scanline * 256 + NST_MIN(cycles.hClock,255) : ~0U;
		}

		NST_FORCE_INLINE bool Ppu::IsDead() const
		{
			return scanline == SCANLINE_VBLANK || !(regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED);
		}

		NST_FORCE_INLINE uint Ppu::Coloring() const
		{
			return (regs.ctrl[1] & Regs::CTRL1_MONOCHROME) ? Palette::MONO : Palette::COLOR;
		}

		NST_FORCE_INLINE uint Ppu::Emphasis() const
		{
			return (regs.ctrl[1] & Regs::CTRL1_EMPHASIS) << 1;
		}

		NES_POKE_D(Ppu,2000)
		{
			Update( cycles.one );

			NST_VERIFY( cpu.GetCycles() >= cycles.reset || !data );

			if (cpu.GetCycles() >= cycles.reset)
			{
				scroll.latch = (scroll.latch & 0x73FF) | (data & 0x03) << 10;
				oam.height = (data >> 2 & 8) + 8;

				io.latch = data;
				data = regs.ctrl[0] ;
				regs.ctrl[0] = io.latch;

				if ((regs.ctrl[0] & regs.status & Regs::CTRL0_NMI) > data)
				{
					const Cycle clock = cpu.GetCycles() + cycles.one;

					if (clock < GetHVIntClock())
						cpu.DoNMI( clock );
				}
			}
		}

		NES_POKE_D(Ppu,2001)
		{
			Update( cycles.one );

			NST_VERIFY( cpu.GetCycles() >= cycles.reset || !data );

			if (cpu.GetCycles() >= cycles.reset)
			{
				if ((regs.ctrl[1] ^ data) & (Regs::CTRL1_BG_ENABLED_NO_CLIP|Regs::CTRL1_SP_ENABLED_NO_CLIP))
				{
					/*original code begin*/
					//tiles.show[0] = (data & Regs::CTRL1_BG_ENABLED) ? 0xFF : 0x00;
					//tiles.show[1] = (data & Regs::CTRL1_BG_ENABLED_NO_CLIP) == Regs::CTRL1_BG_ENABLED_NO_CLIP ? 0xFF : 0x00;

					//oam.show[0] = (data & Regs::CTRL1_SP_ENABLED) ? 0xFF : 0x00;
					//oam.show[1] = (data & Regs::CTRL1_SP_ENABLED_NO_CLIP) == Regs::CTRL1_SP_ENABLED_NO_CLIP ? 0xFF : 0x00;
					/*original code end*/
					/*skin hack begin*/
					tiles.show[0] = (data & Regs::CTRL1_BG_ENABLED) ? 0xFFFF : 0x00;
					tiles.show[1] = (data & Regs::CTRL1_BG_ENABLED_NO_CLIP) == Regs::CTRL1_BG_ENABLED_NO_CLIP ? 0xFFFF : 0x00;

					oam.show[0] = (data & Regs::CTRL1_SP_ENABLED) ? 0xFFFF : 0x00;
					oam.show[1] = (data & Regs::CTRL1_SP_ENABLED_NO_CLIP) == Regs::CTRL1_SP_ENABLED_NO_CLIP ? 0xFFFF : 0x00;
					/*skin hack end*/
					const uint pos = (cycles.hClock - 8) >= (256-16);

					tiles.mask = tiles.show[pos];
					oam.mask = oam.show[pos];

					if ((regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED) && !(data & Regs::CTRL1_BG_SP_ENABLED))
						UpdateAddressLine(scroll.address & 0x3fff);
				}

				io.latch = data;
				data = (regs.ctrl[1] ^ data) & (Regs::CTRL1_EMPHASIS|Regs::CTRL1_MONOCHROME);
				regs.ctrl[1] = io.latch;

				if (data)
				{
					const uint ce[] = { Coloring(), Emphasis() };

					const byte* const NST_RESTRICT map = rgbMap;

					if (!map)
					{
						for (uint i=0; i < Palette::SIZE; ++i)
							output.palette[i] = palette.ram[i] & ce[0] | ce[1];
					}
					else
					{
						for (uint i=0; i < Palette::SIZE; ++i)
							output.palette[i] = map[palette.ram[i] & Palette::COLOR] & ce[0] | ce[1];
					}
				}
			}
		}

		NES_PEEK_A(Ppu,2002)
		{
			Update( cycles.one, address );

			uint status = regs.status & 0xFF;

			regs.status &= (Regs::STATUS_VBLANK^0xFFU);
			scroll.toggle = 0;
			io.latch = (io.latch & Regs::STATUS_LATCH) | status;

			return io.latch;
		}

		NES_PEEK_A(Ppu,2002_RC2C05_01_04)
		{
			return NES_DO_PEEK(2002,address) & 0xC0 | 0x1B;
		}

		NES_PEEK_A(Ppu,2002_RC2C05_02)
		{
			return NES_DO_PEEK(2002,address) & 0xC0 | 0x3D;
		}

		NES_PEEK_A(Ppu,2002_RC2C05_03)
		{
			return NES_DO_PEEK(2002,address) & 0xC0 | 0x1C;
		}

		NES_POKE_D(Ppu,2003)
		{
			Update( cycles.one );

			regs.oam = data;
			io.latch = data;
		}

		NES_POKE_D(Ppu,2004)
		{
			Update( cycles.one );

			NST_ASSERT( regs.oam < Oam::SIZE );
			NST_VERIFY( IsDead() );

			if (IsDead())
			{
				if ((regs.oam & 0x03) == 0x02)
					data &= 0xE3;
			}
			else
			{
				data = 0xFF;
			}

			byte* const NST_RESTRICT value = oam.ram + regs.oam;
			regs.oam = (regs.oam + 1) & 0xFF;
			io.latch = data;
			*value = data;
		}

		NES_PEEK(Ppu,2004)
		{
			NST_ASSERT( regs.oam <= 0xFF );

			if (!(regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED) || cpu.GetCycles() - (cpu.GetFrameCycles() - (341 * 241) * cycles.one) >= (341 * 240) * cycles.one)
			{
				io.latch = oam.ram[regs.oam];
			}
			else
			{
				Update( cycles.one );

				io.latch = oam.latch;
			}

			return io.latch;
		}

		NES_POKE_D(Ppu,2005)
		{
			Update( cycles.one );

			NST_VERIFY( cpu.GetCycles() >= cycles.reset || !data );

			if (cpu.GetCycles() >= cycles.reset)
			{
				io.latch = data;

				if (scroll.toggle ^= 1)
				{
					scroll.latch = (scroll.latch & 0x7FE0) | (data >> 3);
					scroll.xFine = data & 0x7;
				}
				else
				{
					scroll.latch = (scroll.latch & 0x0C1F) | ((data << 2 | data << 12) & 0x73E0);
				}
			}
		}

		NES_POKE_D(Ppu,2006)
		{
			Update( cycles.one );

			NST_VERIFY( cpu.GetCycles() >= cycles.reset || !data );

			if (cpu.GetCycles() >= cycles.reset)
			{
				io.latch = data;

				if (scroll.toggle ^= 1)
				{
					scroll.latch = (scroll.latch & 0x00FF) | (data & 0x3F) << 8;
				}
				else
				{
					scroll.latch = (scroll.latch & 0x7F00) | data;
					scroll.address = scroll.latch;
					if (!(regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED) ||
					    (scanline == SCANLINE_VBLANK)) {
						UpdateAddressLine(scroll.address & 0x3fff);
					}
				}
			}
		}

		NES_POKE_D(Ppu,2007)
		{
			Update( cycles.one * 4 );

			uint address = scroll.address;

			UpdateVramAddress();
			if (!(regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED) ||
			    (scanline == SCANLINE_VBLANK)) {
				UpdateAddressLine(scroll.address & 0x3fff);
			}
			else {
				return;
			}

			io.latch = data;

			if ((address & 0x3F00) == 0x3F00)
			{
				address &= 0x1F;

				const uint final = (!rgbMap ? data : rgbMap[data & Palette::COLOR]) & Coloring() | Emphasis();

				palette.ram[address] = data;
				output.palette[address] = final;

				if (!(address & 0x3))
				{
					palette.ram[address ^ 0x10] = data;
					output.palette[address ^ 0x10] = final;
				}
				
				output.bgColor = palette.ram[0] & uint(Palette::COLOR);
			}
			else
			{
				address &= 0x3FFF;

				if (address >= 0x2000)
					nmt.Poke( address & 0xFFF, data );
				else
					chr.Poke( address, data );
			}
		}

		NES_PEEK_A(Ppu,2007)
		{
			Update( cycles.one, address );

			address = scroll.address & 0x3FFF;
			UpdateVramAddress();
			if (!(regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED) ||
			    (scanline == SCANLINE_VBLANK)) {
				UpdateAddressLine(scroll.address & 0x3fff);
			}

			io.latch = (address & 0x3F00) != 0x3F00 ? io.buffer : palette.ram[address & 0x1F] & Coloring();
			io.buffer = (address >= 0x2000 ? nmt.FetchName( address ) : chr.FetchPattern( address ));

			return io.latch;
		}

		NES_POKE_D(Ppu,2xxx)
		{
			io.latch = data;
		}

		NES_PEEK(Ppu,2xxx)
		{
			return io.latch;
		}

		NES_PEEK(Ppu,3000)
		{
			Update( cycles.one );

			return io.latch;
		}

		NES_POKE_D(Ppu,4014)
		{
			if (cpu.IsOddCycle())
				cpu.StealCycles( cpu.GetClock() );

			Update( cycles.one );
			cpu.StealCycles( cpu.GetClock() );

			NST_ASSERT( regs.oam < 0x100 );

			data <<= 8;

			if ((regs.oam == 0x00 && data < 0x2000) && (!(regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED) || cpu.GetCycles() <= GetHVIntClock() - cpu.GetClock() * 512))
			{
				cpu.StealCycles( cpu.GetClock() * 512 );

				const byte* const NST_RESTRICT cpuRam = cpu.GetRam() + (data & (Cpu::RAM_SIZE-1));
				byte* const NST_RESTRICT oamRam = oam.ram;

				for (uint i=0x00; i < 0x100; i += 0x4)
				{
					oamRam[i+0x0] = cpuRam[i+0x0];
					oamRam[i+0x1] = cpuRam[i+0x1];
					oamRam[i+0x2] = cpuRam[i+0x2] & 0xE3U;
					oamRam[i+0x3] = cpuRam[i+0x3];
				}

				io.latch = oamRam[0xFF];
			}
			else do
			{
				io.latch = cpu.Peek( data++ );
				cpu.StealCycles( cpu.GetClock() );

				Update( cycles.one );
				cpu.StealCycles( cpu.GetClock() );

				NST_VERIFY( IsDead() );

				if (IsDead())
				{
					if ((regs.oam & 0x03) == 0x02)
						io.latch &= 0xE3;
				}
				else
				{
					io.latch = 0xFF;
				}

				byte* const NST_RESTRICT out = oam.ram + regs.oam;
				regs.oam = (regs.oam + 1) & 0xFF;
				*out = io.latch;
			}
			while (data & 0xFF);
		}

		NES_PEEK(Ppu,4014)
		{
			return 0x40;
		}

		NST_FORCE_INLINE void Ppu::Scroll::ClockX()
		{
			if ((address & X_TILE) != X_TILE)
				address++;
			else
				address ^= (X_TILE|NAME_LOW);
		}

		NST_SINGLE_CALL void Ppu::Scroll::ResetX()
		{
			address = (address & ((X_TILE|NAME_LOW) ^ 0x7FFFU)) | (latch & (X_TILE|NAME_LOW));
		}

		NST_SINGLE_CALL void Ppu::Scroll::ClockY()
		{
			if ((address & Y_FINE) != (7U << 12))
			{
				address += (1U << 12);
			}
			else switch (address & Y_TILE)
			{
				default:         address = (address & (Y_FINE ^ 0x7FFFU)) + (1U << 5); break;
				case (29U << 5): address ^= NAME_HIGH;
				case (31U << 5): address &= (Y_FINE|Y_TILE) ^ 0x7FFFU; break;
			}
		}

		NST_SINGLE_CALL void Ppu::PreLoadTiles()
		{
			if( SkinHack::tryFetchBg(io.pattern, tiles.pattern[0], tiles.pattern[1], tiles.attribute, tiles.pixels, output.palette, cycles.hClock) ) {
				return;
			}

			const byte* const NST_RESTRICT src[] =
			{
				tileLut.block[tiles.pattern[0] | (tiles.attribute & 0x3U) << 8],
				tileLut.block[tiles.pattern[1] | (tiles.attribute & 0x3U) << 8]
			};

			NST_ASSERT( tiles.index == 8 );

			/*original code begin*/
			// byte* const NST_RESTRICT dst = tiles.pixels;
			/*original code end*/
			/*skin hack begin*/
			word* const NST_RESTRICT dst = tiles.pixels;
			/*skin hack end*/

			dst[0] = src[0][0];
			dst[1] = src[1][0];
			dst[2] = src[0][1];
			dst[3] = src[1][1];
			dst[4] = src[0][2];
			dst[5] = src[1][2];
			dst[6] = src[0][3];
			dst[7] = src[1][3];
		}

		NST_SINGLE_CALL void Ppu::LoadTiles()
		{
			const byte* const NST_RESTRICT src[] =
			{
				tileLut.block[tiles.pattern[0] | (tiles.attribute & 0x3U) << 8],
				tileLut.block[tiles.pattern[1] | (tiles.attribute & 0x3U) << 8]
			};

			NST_ASSERT( tiles.index == 0 || tiles.index == 8 );

			/*original code begin*/
			// byte* const NST_RESTRICT dst = tiles.pixels + tiles.index;
			/*original code end*/
			/*skin hack begin*/
			word* const NST_RESTRICT dst = tiles.pixels + tiles.index;
			/*skin hack end*/
			tiles.index ^= 8U;

			if( SkinHack::tryFetchBg(io.pattern, tiles.pattern[0], tiles.pattern[1], tiles.attribute, dst, output.palette, cycles.hClock) ) {
				return;
			}

			dst[0] = src[0][0];
			dst[1] = src[1][0];
			dst[2] = src[0][1];
			dst[3] = src[1][1];
			dst[4] = src[0][2];
			dst[5] = src[1][2];
			dst[6] = src[0][3];
			dst[7] = src[1][3];
		}

		NST_FORCE_INLINE void Ppu::EvaluateSpritesEven()
		{
			if (cycles.hClock >= 64)
				oam.latch = oam.ram[oam.address];
		}

		NST_FORCE_INLINE void Ppu::EvaluateSpritesOdd()
		{
			(*this.*oam.phase)();
		}

		void Ppu::EvaluateSpritesPhase0()
		{
		}

		void Ppu::EvaluateSpritesPhase1()
		{
			oam.index++;

			if (scanline - oam.latch >= oam.height)
			{
				if (oam.index != 64)
				{
					oam.address = (oam.index != 2 ? oam.address + 4 : 8);
				}
				else
				{
					oam.address = 0;
					oam.phase = &Ppu::EvaluateSpritesPhase9;
				}
			}
			else
			{
				oam.address++;
				oam.phase = &Ppu::EvaluateSpritesPhase2;
				oam.buffered[0] = oam.latch;
			}
		}

		void Ppu::EvaluateSpritesPhase2()
		{
			oam.address++;
			oam.phase = &Ppu::EvaluateSpritesPhase3;
			oam.buffered[1] = oam.latch;
		}

		void Ppu::EvaluateSpritesPhase3()
		{
			oam.address++;
			oam.phase = &Ppu::EvaluateSpritesPhase4;
			oam.buffered[2] = oam.latch;
		}

		void Ppu::EvaluateSpritesPhase4()
		{
			oam.buffered[3] = oam.latch;
			oam.buffered += 4;

			if (oam.index != 64)
			{
				oam.phase = (oam.buffered != oam.limit ? &Ppu::EvaluateSpritesPhase1 : &Ppu::EvaluateSpritesPhase5);

				if (oam.index != 2)
				{
					oam.address++;

					if (oam.index == 1)
						oam.spriteZeroInLine = true;
				}
				else
				{
					oam.address = 8;
				}
			}
			else
			{
				oam.address = 0;
				oam.phase = &Ppu::EvaluateSpritesPhase9;
			}
		}

		void Ppu::EvaluateSpritesPhase5()
		{
			if (scanline - oam.latch >= oam.height)
			{
				oam.address = ((oam.address + 4) & 0xFC) + ((oam.address + 1) & 0x03);

				if (oam.address <= 5)
				{
					oam.phase = &Ppu::EvaluateSpritesPhase9;
					oam.address &= 0xFC;
				}
			}
			else
			{
				oam.phase = &Ppu::EvaluateSpritesPhase6;
				oam.address = (oam.address + 1) & 0xFF;
				regs.status |= Regs::STATUS_SP_OVERFLOW;
			}
		}

		void Ppu::EvaluateSpritesPhase6()
		{
			oam.phase = &Ppu::EvaluateSpritesPhase7;
			oam.address = (oam.address + 1) & 0xFF;
		}

		void Ppu::EvaluateSpritesPhase7()
		{
			oam.phase = &Ppu::EvaluateSpritesPhase8;
			oam.address = (oam.address + 1) & 0xFF;
		}

		void Ppu::EvaluateSpritesPhase8()
		{
			oam.phase = &Ppu::EvaluateSpritesPhase9;
			oam.address = (oam.address + 1) & 0xFF;

			if ((oam.address & 0x3) == 0x3)
				oam.address++;

			oam.address &= 0xFC;
		}

		void Ppu::EvaluateSpritesPhase9()
		{
			oam.address = (oam.address + 4) & 0xFF;
		}

		NST_FORCE_INLINE uint Ppu::OpenSprite() const
		{
			return (regs.ctrl[0] & (Regs::CTRL0_SP_OFFSET|Regs::CTRL0_SP8X16)) ? 0x1FF0 : 0x0FF0;
		}

		NST_FORCE_INLINE uint Ppu::OpenSprite(const byte* const NST_RESTRICT buffer) const
		{
			uint address;
			const uint comparitor = (uint(scanline) - buffer[0]) ^ ((buffer[2] & uint(Oam::Y_FLIP)) ? 0xF : 0x0);

			if (regs.ctrl[0] & Regs::CTRL0_SP8X16)
			{
				address =
				(
					((buffer[1] & uint(Oam::TILE_LSB)) << 12) |
					((buffer[1] & (Oam::TILE_LSB ^ 0xFFU)) << 4) |
					((comparitor & Oam::RANGE_MSB) << 1)
				);
			}
			else
			{
				address = (regs.ctrl[0] & Regs::CTRL0_SP_OFFSET) << 9 | buffer[1] << 4;
			}

			return address | comparitor & Oam::XFINE;
		}

		NST_FORCE_INLINE void Ppu::LoadSprite(const uint pattern0,const uint pattern1,const byte* const NST_RESTRICT buffer)
		{
			/*original code begin*/
			//if (pattern0 | pattern1)
			//{
				//uint a = (buffer[2] & uint(Oam::X_FLIP)) ? 7 : 0;

				//uint p =
				//(
					//(pattern0 >> 1 & 0x0055) | (pattern1 << 0 & 0x00AA) |
					//(pattern0 << 8 & 0x5500) | (pattern1 << 9 & 0xAA00)
				//);

				//Oam::Output* const NST_RESTRICT entry = oam.visible++;

				//entry->pixels[( a^=6 )] = ( p       ) & 0x3;
				//entry->pixels[( a^=2 )] = ( p >>= 2 ) & 0x3;
				//entry->pixels[( a^=6 )] = ( p >>= 2 ) & 0x3;
				//entry->pixels[( a^=2 )] = ( p >>= 2 ) & 0x3;
				//entry->pixels[( a^=7 )] = ( p >>= 2 ) & 0x3;
				//entry->pixels[( a^=2 )] = ( p >>= 2 ) & 0x3;
				//entry->pixels[( a^=6 )] = ( p >>= 2 ) & 0x3;
				//entry->pixels[( a^=2 )] = ( p >>= 2 );

				//const uint attribute = buffer[2];

				//entry->x       = buffer[3];
				//entry->palette = Palette::SPRITE_OFFSET + ((attribute & Oam::COLOR) << 2);
				//entry->behind  = (attribute & Oam::BEHIND) ? 0x3 : 0x0;
				//entry->zero    = (buffer == oam.buffer && oam.spriteZeroInLine) ? 0x3 : 0x0;
			//} 			
			/*original code end*/
			/*skin hack begin*/
			uint a = (buffer[2] & uint(Oam::X_FLIP)) ? 7 : 0;

			uint p =
			(
				(pattern0 >> 1 & 0x0055) | (pattern1 << 0 & 0x00AA) |
				(pattern0 << 8 & 0x5500) | (pattern1 << 9 & 0xAA00)
			);

			Oam::Output* const NST_RESTRICT entry = oam.visible;
			if( !SkinHack::tryFetchSprite(io.address, p, buffer[2] & Oam::COLOR, buffer[2] & Oam::X_FLIP, entry->pixels, output.palette) ) {
				if (pattern0 | pattern1) {
					entry->pixels[( a^=6 )] = ( p       ) & 0x3;
					entry->pixels[( a^=2 )] = ( p >>= 2 ) & 0x3;
					entry->pixels[( a^=6 )] = ( p >>= 2 ) & 0x3;
					entry->pixels[( a^=2 )] = ( p >>= 2 ) & 0x3;
					entry->pixels[( a^=7 )] = ( p >>= 2 ) & 0x3;
					entry->pixels[( a^=2 )] = ( p >>= 2 ) & 0x3;
					entry->pixels[( a^=6 )] = ( p >>= 2 ) & 0x3;
					entry->pixels[( a^=2 )] = ( p >>= 2 );
				} else {
					return;
				}
			}
			oam.visible++;
			const uint attribute = buffer[2];

			entry->x       = buffer[3];
			entry->palette = Palette::SPRITE_OFFSET + ((attribute & Oam::COLOR) << 2);
			entry->behind  = (attribute & Oam::BEHIND) ? 0x3 : 0x0;
			entry->zero    = (buffer == oam.buffer && oam.spriteZeroInLine) ? 0x3 : 0x0;
			/*skin hack end*/
		}

		void Ppu::LoadExtendedSprites()
		{
			const byte* NST_RESTRICT buffer = oam.buffer + (8*4);
			NST_ASSERT( buffer < oam.buffered );

			do
			{
				const uint address = OpenSprite( buffer );

				const uint patterns[2] =
				{
					chr.FetchPattern( address | 0x0 ),
					chr.FetchPattern( address | 0x8 )
				};
				/*skin hack begin*/
				uint tmp = io.address;
				io.address = address;
				/*skin hack end*/
				LoadSprite( patterns[0], patterns[1], buffer );
				/*skin hack begin*/
				io.address = tmp;
				/*skin hack end*/
				buffer += 4;
			}
			while (buffer != oam.buffered);
		}

		NST_FORCE_INLINE void Ppu::RenderPixel()
		{
			uint clock;
			uint pixel = tiles.pixels[((clock=cycles.hClock++) + scroll.xFine) & 15] & tiles.mask;

			for (const Oam::Output* NST_RESTRICT sprite=oam.output, *const end=oam.visible; sprite != end; ++sprite)
			{
				uint x = clock - sprite->x;

				if (x > 7)
					continue;

				x = sprite->pixels[x] & oam.mask;

				if (x)
				{
					/*original code begin*/
					//if (pixel & sprite->zero)
					/*original code end*/
					/*skin hack begin*/
					if (pixel && sprite->zero)
					/*skin hack end*/
						regs.status |= Regs::STATUS_SP_ZERO_HIT;

					/*original code begin*/
					//if (!(pixel & sprite->behind))
					//	pixel = sprite->palette + x;
					/*original code end*/
					/*skin hack begin*/
					if (!(pixel && sprite->behind))
						pixel = x & 0x8000 ? x : sprite->palette + x;
					/*skin hack end*/

					break;
				}
			}

			Video::Screen::Pixel* const NST_RESTRICT target = output.target++;
			/*original code begin*/
			// *target = output.palette[pixel];
			/*original code end*/
			/*skin hack begin*/
			*target = pixel & 0x8000 ? pixel : output.palette[pixel];
			/*skin hack end*/
		}

		NST_SINGLE_CALL void Ppu::RenderPixel255()
		{
			cycles.hClock = 256;
			uint pixel = tiles.pixels[(255 + scroll.xFine) & 15] & tiles.mask;

			for (const Oam::Output* NST_RESTRICT sprite=oam.output, *const end=oam.visible; sprite != end; ++sprite)
			{
				uint x = 255U - sprite->x;

				if (x > 7)
					continue;

				x = sprite->pixels[x] & oam.mask;

				if (x)
				{
					/*original code begin*/
					//if (!(pixel & sprite->behind))
					//	pixel = sprite->palette + x;
					/*original code end*/
					/*skin hack begin*/
					if (!(pixel && sprite->behind))
						pixel = x & 0x8000 ? x : sprite->palette + x;
					/*skin hack end*/

					break;
				}
			}

			Video::Screen::Pixel* const NST_RESTRICT target = output.target++;
			/*original code begin*/
			// *target = output.palette[pixel];
			/*original code end*/
			/*skin hack begin*/
			*target = pixel & 0x8000 ? pixel : output.palette[pixel];
			/*skin hack end*/
		}

		NST_NO_INLINE void Ppu::Run()
		{
			NST_VERIFY( cycles.count != cycles.hClock );

			if (regs.ctrl[1] & Regs::CTRL1_BG_SP_ENABLED)
			{
				switch (cycles.hClock)
				{
					case 0:
					case 8:
					case 16:
					case 24:
					case 32:
					case 40:
					case 48:
					case 56:
					case 72:
					case 80:
					case 88:
					case 96:
					case 104:
					case 112:
					case 120:
					case 128:
					case 136:
					case 144:
					case 152:
					case 160:
					case 168:
					case 176:
					case 184:
					case 192:
					case 200:
					case 208:
					case 216:
					case 224:
					case 232:
					case 240:
					case 248:
					HActive:

						LoadTiles();
						EvaluateSpritesEven();
						OpenName();
						RenderPixel();

						if (cycles.count <= cycles.hClock)
							break;

					case 1:
					case 9:
					case 17:
					case 25:
					case 33:
					case 41:
					case 49:
					case 57:
					case 65:
					case 73:
					case 81:
					case 89:
					case 97:
					case 105:
					case 113:
					case 121:
					case 129:
					case 137:
					case 145:
					case 153:
					case 161:
					case 169:
					case 177:
					case 185:
					case 193:
					case 201:
					case 209:
					case 217:
					case 225:
					case 233:
					case 241:
					case 249:

						FetchName();
						EvaluateSpritesOdd();
						RenderPixel();

						if (cycles.count <= cycles.hClock)
							break;

					case 2:
					case 10:
					case 18:
					case 26:
					case 34:
					case 42:
					case 50:
					case 58:
					case 66:
					case 74:
					case 82:
					case 90:
					case 98:
					case 106:
					case 114:
					case 122:
					case 130:
					case 138:
					case 146:
					case 154:
					case 162:
					case 170:
					case 178:
					case 186:
					case 194:
					case 202:
					case 210:
					case 218:
					case 226:
					case 234:
					case 242:
					case 250:

						EvaluateSpritesEven();
						OpenAttribute();
						RenderPixel();

						if (cycles.count <= cycles.hClock)
							break;

					case 3:
					case 11:
					case 19:
					case 27:
					case 35:
					case 43:
					case 51:
					case 59:
					case 67:
					case 75:
					case 83:
					case 91:
					case 99:
					case 107:
					case 115:
					case 123:
					case 131:
					case 139:
					case 147:
					case 155:
					case 163:
					case 171:
					case 179:
					case 187:
					case 195:
					case 203:
					case 211:
					case 219:
					case 227:
					case 235:
					case 243:
					case 251:

						FetchAttribute();
						EvaluateSpritesOdd();

						if (cycles.hClock == 251)
							scroll.ClockY();

						scroll.ClockX();
						RenderPixel();

						if (cycles.count <= cycles.hClock)
							break;

					case 4:
					case 12:
					case 20:
					case 28:
					case 36:
					case 44:
					case 52:
					case 60:
					case 68:
					case 76:
					case 84:
					case 92:
					case 100:
					case 108:
					case 116:
					case 124:
					case 132:
					case 140:
					case 148:
					case 156:
					case 164:
					case 172:
					case 180:
					case 188:
					case 196:
					case 204:
					case 212:
					case 220:
					case 228:
					case 236:
					case 244:
					case 252:

						EvaluateSpritesEven();
						OpenPattern( io.pattern | 0x0 );
						RenderPixel();

						if (cycles.count <= cycles.hClock)
							break;

					case 5:
					case 13:
					case 21:
					case 29:
					case 37:
					case 45:
					case 53:
					case 61:
					case 69:
					case 77:
					case 85:
					case 93:
					case 101:
					case 109:
					case 117:
					case 125:
					case 133:
					case 141:
					case 149:
					case 157:
					case 165:
					case 173:
					case 181:
					case 189:
					case 197:
					case 205:
					case 213:
					case 221:
					case 229:
					case 237:
					case 245:
					case 253:

						FetchBgPattern0();
						EvaluateSpritesOdd();
						RenderPixel();

						if (cycles.count <= cycles.hClock)
							break;

					case 6:
					case 14:
					case 22:
					case 30:
					case 38:
					case 46:
					case 54:
					case 62:
					case 70:
					case 78:
					case 86:
					case 94:
					case 102:
					case 110:
					case 118:
					case 126:
					case 134:
					case 142:
					case 150:
					case 158:
					case 166:
					case 174:
					case 182:
					case 190:
					case 198:
					case 206:
					case 214:
					case 222:
					case 230:
					case 238:
					case 246:
					case 254:

						EvaluateSpritesEven();
						OpenPattern( io.pattern | 0x8 );
						RenderPixel();

						if (cycles.count <= cycles.hClock)
							break;

						if (cycles.hClock == 255)
							goto HActive255;

					case 7:
					case 15:
					case 23:
					case 31:
					case 39:
					case 47:
					case 55:
					case 63:
					case 71:
					case 79:
					case 87:
					case 95:
					case 103:
					case 111:
					case 119:
					case 127:
					case 135:
					case 143:
					case 151:
					case 159:
					case 167:
					case 175:
					case 183:
					case 191:
					case 199:
					case 207:
					case 215:
					case 223:
					case 231:
					case 239:
					case 247:

						FetchBgPattern1();
						EvaluateSpritesOdd();
						RenderPixel();
						tiles.mask = tiles.show[0];
						oam.mask = oam.show[0];

						if (cycles.count <= cycles.hClock)
							break;

						if (cycles.hClock != 64)
							goto HActive;

					case 64:

						NST_VERIFY( regs.oam == 0 );
						oam.address = regs.oam & Oam::OFFSET_TO_0_1;
						oam.phase = &Ppu::EvaluateSpritesPhase1;
						oam.latch = 0xFF;
						goto HActive;

					case 255:
					HActive255:

						FetchBgPattern1();
						EvaluateSpritesOdd();
						RenderPixel255();

						if (cycles.count <= 256)
							break;

					case 256:

						OpenName();
						oam.latch = 0xFF;
						cycles.hClock = 257;

						if (cycles.count <= 257)
							break;

					case 257:

						if (hBlankHook)
							hBlankHook.Execute();

						scroll.ResetX();
						oam.visible = oam.output;
						cycles.hClock = 258;

						if (cycles.count <= 258)
							break;

					case 258:
					case 266:
					case 274:
					case 282:
					case 290:
					case 298:
					case 306:
					case 314:
					HBlankSp:

						OpenName();
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

					case 260:
					case 268:
					case 276:
					case 284:
					case 292:
					case 300:
					case 308:
					case 316:
					{
						const byte* const buffer = oam.buffer + ((cycles.hClock - 260) / 2);
						OpenPattern( buffer >= oam.buffered ? OpenSprite() : OpenSprite(buffer) );

						if (scanline == 238 && cycles.hClock == 316)
							regs.oam = 0;

						if (cycles.count <= ++cycles.hClock)
							break;
					}

					case 261:
					case 269:
					case 277:
					case 285:
					case 293:
					case 301:
					case 309:
					case 317:

						if (oam.buffer + ((cycles.hClock - 261) / 2) < oam.buffered)
							io.pattern = FetchSpPattern();

						if (cycles.count <= ++cycles.hClock)
							break;

					case 262:
					case 270:
					case 278:
					case 286:
					case 294:
					case 302:
					case 310:
					case 318:

						OpenPattern( io.address | 0x8 );

						if (cycles.count <= ++cycles.hClock)
							break;

					case 263:
					case 271:
					case 279:
					case 287:
					case 295:
					case 303:
					case 311:
					case 319:
					{
						const byte* const buffer = oam.buffer + ((cycles.hClock - 263) / 2);

						if (buffer < oam.buffered)
							LoadSprite( io.pattern, FetchSpPattern(), buffer );

						if (cycles.count <= ++cycles.hClock)
							break;

						if (cycles.hClock == 320)
							goto HBlankBg;
					}

					case 264:
					case 272:
					case 280:
					case 288:
					case 296:
					case 304:
					case 312:

						OpenName();
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

						goto HBlankSp;

					case 320:
					HBlankBg:

						if (oam.buffer + (8*4) < oam.buffered)
							LoadExtendedSprites();

						OpenName();

						if (hActiveHook)
							hActiveHook.Execute();

						oam.latch = oam.ram[0];
						oam.buffered = oam.buffer;
						oam.spriteZeroInLine = false;
						oam.index = 0;
						oam.phase = &Ppu::EvaluateSpritesPhase0;
						cycles.hClock = 321;

						if (cycles.count <= 321)
							break;

					case 321:

						FetchName();
						cycles.hClock = 322;

						if (cycles.count <= 322)
							break;

					case 322:

						OpenAttribute();
						cycles.hClock = 323;

						if (cycles.count <= 323)
							break;

					case 323:

						FetchAttribute();
						scroll.ClockX();
						cycles.hClock = 324;

						if (cycles.count <= 324)
							break;

					case 324:

						OpenPattern( io.pattern | 0x0 );
						cycles.hClock = 325;

						if (cycles.count <= 325)
							break;

					case 325:

						FetchBgPattern0();
						cycles.hClock = 326;

						if (cycles.count <= 326)
							break;

					case 326:

						OpenPattern( io.pattern | 0x8 );
						cycles.hClock = 327;

						if (cycles.count <= 327)
							break;

					case 327:

						FetchBgPattern1();
						cycles.hClock = 328;

						if (cycles.count <= 328)
							break;

					case 328:

						PreLoadTiles();
						OpenName();
						cycles.hClock = 329;

						if (cycles.count <= 329)
							break;

					case 329:

						FetchName();
						cycles.hClock = 330;

						if (cycles.count <= 330)
							break;

					case 330:

						OpenAttribute();
						cycles.hClock = 331;

						if (cycles.count <= 331)
							break;

					case 331:

						FetchAttribute();
						scroll.ClockX();
						cycles.hClock = 332;

						if (cycles.count <= 332)
							break;

					case 332:

						OpenPattern( io.pattern | 0x0 );
						cycles.hClock = 333;

						if (cycles.count <= 333)
							break;

					case 333:

						FetchBgPattern0();
						cycles.hClock = 334;

						if (cycles.count <= 334)
							break;

					case 334:

						OpenPattern( io.pattern | 0x8 );
						cycles.hClock = 335;

						if (cycles.count <= 335)
							break;

					case 335:

						FetchBgPattern1();
						cycles.hClock = 336;

						if (cycles.count <= 336)
							break;

					case 336:

						OpenName();
						cycles.hClock = 337;

						if (cycles.count <= 337)
							break;

					case 337:

						tiles.mask = tiles.show[1];
						oam.mask = oam.show[1];

						if (scanline == SCANLINE_HDUMMY && model == PPU_RP2C02)
						{
							if (regs.frame)
							{
								output.burstPhase = (output.burstPhase + 2) % 3;
								cpu.SetFrameCycles( PPU_RP2C02_HVSYNC_1 );
							}
							else
							{
								output.burstPhase = (output.burstPhase + 1) % 3;
							}
						}

						cycles.hClock = 338;

						if (cycles.count <= 338)
							break;

					case 338:

						OpenName();

						if (scanline++ != 239)
						{
							const uint line = (scanline != 0 || model != PPU_RP2C02 || !regs.frame ? 341 : 340);

							cycles.hClock = 0;
							cycles.vClock += line;

							if (cycles.count <= line)
								break;

							cycles.count -= line;

							goto HActive;
						}
						else
						{
							cycles.hClock = HCLOCK_POSTRENDER;

							if (cycles.count <= HCLOCK_POSTRENDER)
								break;
						}

					case HCLOCK_POSTRENDER:
						UpdateAddressLine(scroll.address & 0x3fff);
						cycles.hClock = HCLOCK_VBLANK_0;

						if (cycles.count <= HCLOCK_VBLANK_0)
							break;

					case HCLOCK_VBLANK_0:
					VBlank0:

						regs.status |= Regs::STATUS_VBLANKING;
						cycles.hClock = HCLOCK_VBLANK_1;

						if (cycles.count <= HCLOCK_VBLANK_1)
							break;

					case HCLOCK_VBLANK_1:
					VBlank1:

						regs.status = (regs.status & 0xFF) | (regs.status >> 1 & Regs::STATUS_VBLANK);
						oam.visible = oam.output;
						cycles.hClock = HCLOCK_VBLANK_2;

						if (cycles.count <= HCLOCK_VBLANK_2)
							break;

					case HCLOCK_VBLANK_2:
					VBlank2:

						cycles.hClock = HCLOCK_DUMMY;
						cycles.count = Cpu::CYCLE_MAX;
						cycles.reset = 0;

						if (regs.ctrl[0] & regs.status & Regs::CTRL0_NMI)
							cpu.DoNMI( cpu.GetFrameCycles() );

						return;

					case HCLOCK_BOOT:
						goto Boot;

					case HCLOCK_DUMMY+0:

						regs.status = 0;
						scanline = SCANLINE_HDUMMY;

					case HCLOCK_DUMMY+8:
					case HCLOCK_DUMMY+16:
					case HCLOCK_DUMMY+24:
					case HCLOCK_DUMMY+32:
					case HCLOCK_DUMMY+40:
					case HCLOCK_DUMMY+48:
					case HCLOCK_DUMMY+56:
					case HCLOCK_DUMMY+64:
					case HCLOCK_DUMMY+72:
					case HCLOCK_DUMMY+80:
					case HCLOCK_DUMMY+88:
					case HCLOCK_DUMMY+96:
					case HCLOCK_DUMMY+104:
					case HCLOCK_DUMMY+112:
					case HCLOCK_DUMMY+120:
					case HCLOCK_DUMMY+128:
					case HCLOCK_DUMMY+136:
					case HCLOCK_DUMMY+144:
					case HCLOCK_DUMMY+152:
					case HCLOCK_DUMMY+160:
					case HCLOCK_DUMMY+168:
					case HCLOCK_DUMMY+176:
					case HCLOCK_DUMMY+184:
					case HCLOCK_DUMMY+192:
					case HCLOCK_DUMMY+200:
					case HCLOCK_DUMMY+208:
					case HCLOCK_DUMMY+216:
					case HCLOCK_DUMMY+224:
					case HCLOCK_DUMMY+232:
					case HCLOCK_DUMMY+240:
					case HCLOCK_DUMMY+248:
					HDummyBg:

						OpenName();
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

					case HCLOCK_DUMMY+2:
					case HCLOCK_DUMMY+10:
					case HCLOCK_DUMMY+18:
					case HCLOCK_DUMMY+26:
					case HCLOCK_DUMMY+34:
					case HCLOCK_DUMMY+42:
					case HCLOCK_DUMMY+50:
					case HCLOCK_DUMMY+58:
					case HCLOCK_DUMMY+66:
					case HCLOCK_DUMMY+74:
					case HCLOCK_DUMMY+82:
					case HCLOCK_DUMMY+90:
					case HCLOCK_DUMMY+98:
					case HCLOCK_DUMMY+106:
					case HCLOCK_DUMMY+114:
					case HCLOCK_DUMMY+122:
					case HCLOCK_DUMMY+130:
					case HCLOCK_DUMMY+138:
					case HCLOCK_DUMMY+146:
					case HCLOCK_DUMMY+154:
					case HCLOCK_DUMMY+162:
					case HCLOCK_DUMMY+170:
					case HCLOCK_DUMMY+178:
					case HCLOCK_DUMMY+186:
					case HCLOCK_DUMMY+194:
					case HCLOCK_DUMMY+202:
					case HCLOCK_DUMMY+210:
					case HCLOCK_DUMMY+218:
					case HCLOCK_DUMMY+226:
					case HCLOCK_DUMMY+234:
					case HCLOCK_DUMMY+242:
					case HCLOCK_DUMMY+250:

						OpenAttribute();
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

					case HCLOCK_DUMMY+4:
					case HCLOCK_DUMMY+12:
					case HCLOCK_DUMMY+20:
					case HCLOCK_DUMMY+28:
					case HCLOCK_DUMMY+36:
					case HCLOCK_DUMMY+44:
					case HCLOCK_DUMMY+52:
					case HCLOCK_DUMMY+60:
					case HCLOCK_DUMMY+68:
					case HCLOCK_DUMMY+76:
					case HCLOCK_DUMMY+84:
					case HCLOCK_DUMMY+92:
					case HCLOCK_DUMMY+100:
					case HCLOCK_DUMMY+108:
					case HCLOCK_DUMMY+116:
					case HCLOCK_DUMMY+124:
					case HCLOCK_DUMMY+132:
					case HCLOCK_DUMMY+140:
					case HCLOCK_DUMMY+148:
					case HCLOCK_DUMMY+156:
					case HCLOCK_DUMMY+164:
					case HCLOCK_DUMMY+172:
					case HCLOCK_DUMMY+180:
					case HCLOCK_DUMMY+188:
					case HCLOCK_DUMMY+196:
					case HCLOCK_DUMMY+204:
					case HCLOCK_DUMMY+212:
					case HCLOCK_DUMMY+220:
					case HCLOCK_DUMMY+228:
					case HCLOCK_DUMMY+236:
					case HCLOCK_DUMMY+244:
					case HCLOCK_DUMMY+252:

						OpenPattern( regs.ctrl[0] << 8 & 0x1000 );
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

					case HCLOCK_DUMMY+6:
					case HCLOCK_DUMMY+14:
					case HCLOCK_DUMMY+22:
					case HCLOCK_DUMMY+30:
					case HCLOCK_DUMMY+38:
					case HCLOCK_DUMMY+46:
					case HCLOCK_DUMMY+54:
					case HCLOCK_DUMMY+62:
					case HCLOCK_DUMMY+70:
					case HCLOCK_DUMMY+78:
					case HCLOCK_DUMMY+86:
					case HCLOCK_DUMMY+94:
					case HCLOCK_DUMMY+102:
					case HCLOCK_DUMMY+110:
					case HCLOCK_DUMMY+118:
					case HCLOCK_DUMMY+126:
					case HCLOCK_DUMMY+134:
					case HCLOCK_DUMMY+142:
					case HCLOCK_DUMMY+150:
					case HCLOCK_DUMMY+158:
					case HCLOCK_DUMMY+166:
					case HCLOCK_DUMMY+174:
					case HCLOCK_DUMMY+182:
					case HCLOCK_DUMMY+190:
					case HCLOCK_DUMMY+198:
					case HCLOCK_DUMMY+206:
					case HCLOCK_DUMMY+214:
					case HCLOCK_DUMMY+222:
					case HCLOCK_DUMMY+230:
					case HCLOCK_DUMMY+238:
					case HCLOCK_DUMMY+246:
					case HCLOCK_DUMMY+254:

						OpenPattern( io.address | 0x8 );
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

						if (cycles.hClock != HCLOCK_DUMMY+256)
							goto HDummyBg;

					case HCLOCK_DUMMY+256:
					case HCLOCK_DUMMY+264:
					case HCLOCK_DUMMY+272:
					case HCLOCK_DUMMY+280:
					case HCLOCK_DUMMY+288:
					case HCLOCK_DUMMY+296:
					case HCLOCK_DUMMY+312:
					HDummySp:

						OpenName();
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

					case HCLOCK_DUMMY+258:
					case HCLOCK_DUMMY+266:
					case HCLOCK_DUMMY+274:
					case HCLOCK_DUMMY+282:
					case HCLOCK_DUMMY+290:
					case HCLOCK_DUMMY+298:
					case HCLOCK_DUMMY+306:
					case HCLOCK_DUMMY+314:

						OpenName();
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

					case HCLOCK_DUMMY+260:
					case HCLOCK_DUMMY+268:
					case HCLOCK_DUMMY+276:
					case HCLOCK_DUMMY+284:
					case HCLOCK_DUMMY+292:
					case HCLOCK_DUMMY+300:
					case HCLOCK_DUMMY+308:
					case HCLOCK_DUMMY+316:

						OpenPattern( OpenSprite() );
						cycles.hClock += 2;

						if (cycles.count <= cycles.hClock)
							break;

					case HCLOCK_DUMMY+262:
					case HCLOCK_DUMMY+270:
					case HCLOCK_DUMMY+278:
					case HCLOCK_DUMMY+286:
					case HCLOCK_DUMMY+294:
					case HCLOCK_DUMMY+302:
					case HCLOCK_DUMMY+310:
					case HCLOCK_DUMMY+318:

						OpenPattern( io.address | 0x8 );

						if (cycles.hClock != HCLOCK_DUMMY+318)
						{
							cycles.hClock += 2;

							if (cycles.count <= cycles.hClock)
								break;

							if (cycles.hClock != HCLOCK_DUMMY+304)
								goto HDummySp;
						}
						else
						{
							cycles.hClock = 320;
							cycles.vClock += HCLOCK_DUMMY;
							cycles.count -= HCLOCK_DUMMY;

							if (cycles.count <= cycles.hClock)
								break;

							goto HBlankBg;
						}

					case HCLOCK_DUMMY+304:

						scroll.address = scroll.latch;
						goto HDummySp;

					default:

						NST_UNREACHABLE();
				}
			}
			else
			{
				switch (cycles.hClock)
				{
					case 0:
					case 1:
					case 2:
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
					case 9:
					case 10:
					case 11:
					case 12:
					case 13:
					case 14:
					case 15:
					case 16:
					case 17:
					case 18:
					case 19:
					case 20:
					case 21:
					case 22:
					case 23:
					case 24:
					case 25:
					case 26:
					case 27:
					case 28:
					case 29:
					case 30:
					case 31:
					case 32:
					case 33:
					case 34:
					case 35:
					case 36:
					case 37:
					case 38:
					case 39:
					case 40:
					case 41:
					case 42:
					case 43:
					case 44:
					case 45:
					case 46:
					case 47:
					case 48:
					case 49:
					case 50:
					case 51:
					case 52:
					case 53:
					case 54:
					case 55:
					case 56:
					case 57:
					case 58:
					case 59:
					case 60:
					case 61:
					case 62:
					case 63:
					case 64:
					case 65:
					case 66:
					case 67:
					case 68:
					case 69:
					case 70:
					case 71:
					case 72:
					case 73:
					case 74:
					case 75:
					case 76:
					case 77:
					case 78:
					case 79:
					case 80:
					case 81:
					case 82:
					case 83:
					case 84:
					case 85:
					case 86:
					case 87:
					case 88:
					case 89:
					case 90:
					case 91:
					case 92:
					case 93:
					case 94:
					case 95:
					case 96:
					case 97:
					case 98:
					case 99:
					case 100:
					case 101:
					case 102:
					case 103:
					case 104:
					case 105:
					case 106:
					case 107:
					case 108:
					case 109:
					case 110:
					case 111:
					case 112:
					case 113:
					case 114:
					case 115:
					case 116:
					case 117:
					case 118:
					case 119:
					case 120:
					case 121:
					case 122:
					case 123:
					case 124:
					case 125:
					case 126:
					case 127:
					case 128:
					case 129:
					case 130:
					case 131:
					case 132:
					case 133:
					case 134:
					case 135:
					case 136:
					case 137:
					case 138:
					case 139:
					case 140:
					case 141:
					case 142:
					case 143:
					case 144:
					case 145:
					case 146:
					case 147:
					case 148:
					case 149:
					case 150:
					case 151:
					case 152:
					case 153:
					case 154:
					case 155:
					case 156:
					case 157:
					case 158:
					case 159:
					case 160:
					case 161:
					case 162:
					case 163:
					case 164:
					case 165:
					case 166:
					case 167:
					case 168:
					case 169:
					case 170:
					case 171:
					case 172:
					case 173:
					case 174:
					case 175:
					case 176:
					case 177:
					case 178:
					case 179:
					case 180:
					case 181:
					case 182:
					case 183:
					case 184:
					case 185:
					case 186:
					case 187:
					case 188:
					case 189:
					case 190:
					case 191:
					case 192:
					case 193:
					case 194:
					case 195:
					case 196:
					case 197:
					case 198:
					case 199:
					case 200:
					case 201:
					case 202:
					case 203:
					case 204:
					case 205:
					case 206:
					case 207:
					case 208:
					case 209:
					case 210:
					case 211:
					case 212:
					case 213:
					case 214:
					case 215:
					case 216:
					case 217:
					case 218:
					case 219:
					case 220:
					case 221:
					case 222:
					case 223:
					case 224:
					case 225:
					case 226:
					case 227:
					case 228:
					case 229:
					case 230:
					case 231:
					case 232:
					case 233:
					case 234:
					case 235:
					case 236:
					case 237:
					case 238:
					case 239:
					case 240:
					case 241:
					case 242:
					case 243:
					case 244:
					case 245:
					case 246:
					case 247:
					case 248:
					case 249:
					case 250:
					case 251:
					case 252:
					case 253:
					case 254:
					case 255:
					HActiveOff:
					{
						const uint pixel = output.palette[(scroll.address & 0x3F00) == 0x3F00 ? (scroll.address & 0x001F) : 0];

						uint i = cycles.hClock;
						const uint hClock = NST_MIN(cycles.count,256);
						NST_ASSERT( i < hClock );

						cycles.hClock = hClock;
						tiles.index = (hClock - 1) & 8;

						/*original code begin*/
						// byte* const NST_RESTRICT tile = tiles.pixels;
						/*original code end*/
						/*skin hack begin*/
						word* const NST_RESTRICT tile = tiles.pixels;
						/*skin hack end*/
						Video::Screen::Pixel* NST_RESTRICT target = output.target;

						do
						{
							tile[i++ & 15] = 0;
							*target++ = pixel;
						}
						while (i != hClock);

						output.target = target;

						if (cycles.count <= 256)
							break;
					}

					case 256:

						cycles.hClock = 257;

						if (cycles.count <= 257)
							break;

					case 257:

						if (hBlankHook)
							hBlankHook.Execute();

						oam.visible = oam.output;
						cycles.hClock = 258;

						if (cycles.count <= 258)
							break;

					case 258:
					case 260:
					case 261:
					case 262:
					case 263:
					case 264:
					case 266:
					case 268:
					case 269:
					case 270:
					case 271:
					case 272:
					case 274:
					case 276:
					case 277:
					case 278:
					case 279:
					case 280:
					case 282:
					case 284:
					case 285:
					case 286:
					case 287:
					case 288:
					case 290:
					case 292:
					case 293:
					case 294:
					case 295:
					case 296:
					case 298:
					case 300:
					case 301:
					case 302:
					case 303:
					case 304:
					case 306:
					case 308:
					case 309:
					case 310:
					case 311:
					case 312:
					case 314:
					case 316:
					case 317:
					case 318:
					case 319:

						if (cycles.count <= 320)
						{
							cycles.hClock = cycles.count + ((cycles.count & 0x7) == 3 || (cycles.count & 0x7) == 1);
							break;
						}

					case 320:
					HBlankOff:

						cycles.hClock = 321;

						if (hActiveHook)
							hActiveHook.Execute();

						oam.buffered = oam.buffer;
						oam.spriteZeroInLine = false;
						oam.index = 0;
						oam.phase = &Ppu::EvaluateSpritesPhase0;

						if (cycles.count <= 321)
							break;

					case 321:
					case 322:
					case 323:
					case 324:
					case 325:
					case 326:
					case 327:
					case 328:
					case 329:
					case 330:
					case 331:
					case 332:
					case 333:
					case 334:
					case 335:
					case 336:
					case 337:

						cycles.hClock = cycles.count;

						if (cycles.count <= 338)
							break;

					case 338:

						if (scanline++ != 239)
						{
							tiles.mask = tiles.show[1];
							oam.mask = oam.show[1];

							if (scanline == 0 && model == PPU_RP2C02)
								output.burstPhase = (output.burstPhase + 1) % 3;

							cycles.vClock += 341;
							cycles.hClock = 0;

							if (cycles.count <= 341)
								break;

							cycles.count -= 341;

							goto HActiveOff;
						}
						else
						{
							cycles.hClock = HCLOCK_VBLANK_0;

							if (cycles.count <= HCLOCK_VBLANK_0)
								break;
						}

					case HCLOCK_VBLANK_0:
						goto VBlank0;

					case HCLOCK_VBLANK_1:
						goto VBlank1;

					case HCLOCK_VBLANK_2:
						goto VBlank2;

					case HCLOCK_BOOT:
					Boot:

						regs.status |= Regs::STATUS_VBLANK;
						cycles.hClock = HCLOCK_DUMMY;
						cycles.count = Cpu::CYCLE_MAX;

						if (cycles.reset)
						{
							switch (model)
							{
								case PPU_RP2C07: cycles.reset = PPU_RP2C07_HVREGBOOT - PPU_RP2C07_HVSYNCBOOT; break;
								case PPU_DENDY:  cycles.reset = PPU_DENDY_HVREGBOOT  - PPU_DENDY_HVSYNCBOOT;  break;
								default:         cycles.reset = PPU_RP2C02_HVREGBOOT - PPU_RP2C02_HVSYNCBOOT; break;
							}
						}
						return;

					case HCLOCK_DUMMY+0:

						regs.status = 0;
						scanline = SCANLINE_HDUMMY;

					case HCLOCK_DUMMY+2:
					case HCLOCK_DUMMY+4:
					case HCLOCK_DUMMY+6:
					case HCLOCK_DUMMY+8:
					case HCLOCK_DUMMY+10:
					case HCLOCK_DUMMY+12:
					case HCLOCK_DUMMY+14:
					case HCLOCK_DUMMY+16:
					case HCLOCK_DUMMY+18:
					case HCLOCK_DUMMY+20:
					case HCLOCK_DUMMY+22:
					case HCLOCK_DUMMY+24:
					case HCLOCK_DUMMY+26:
					case HCLOCK_DUMMY+28:
					case HCLOCK_DUMMY+30:
					case HCLOCK_DUMMY+32:
					case HCLOCK_DUMMY+34:
					case HCLOCK_DUMMY+36:
					case HCLOCK_DUMMY+38:
					case HCLOCK_DUMMY+40:
					case HCLOCK_DUMMY+42:
					case HCLOCK_DUMMY+44:
					case HCLOCK_DUMMY+46:
					case HCLOCK_DUMMY+48:
					case HCLOCK_DUMMY+50:
					case HCLOCK_DUMMY+52:
					case HCLOCK_DUMMY+54:
					case HCLOCK_DUMMY+56:
					case HCLOCK_DUMMY+58:
					case HCLOCK_DUMMY+60:
					case HCLOCK_DUMMY+62:
					case HCLOCK_DUMMY+64:
					case HCLOCK_DUMMY+66:
					case HCLOCK_DUMMY+68:
					case HCLOCK_DUMMY+70:
					case HCLOCK_DUMMY+72:
					case HCLOCK_DUMMY+74:
					case HCLOCK_DUMMY+76:
					case HCLOCK_DUMMY+78:
					case HCLOCK_DUMMY+80:
					case HCLOCK_DUMMY+82:
					case HCLOCK_DUMMY+84:
					case HCLOCK_DUMMY+86:
					case HCLOCK_DUMMY+88:
					case HCLOCK_DUMMY+90:
					case HCLOCK_DUMMY+92:
					case HCLOCK_DUMMY+94:
					case HCLOCK_DUMMY+96:
					case HCLOCK_DUMMY+98:
					case HCLOCK_DUMMY+100:
					case HCLOCK_DUMMY+102:
					case HCLOCK_DUMMY+104:
					case HCLOCK_DUMMY+106:
					case HCLOCK_DUMMY+108:
					case HCLOCK_DUMMY+110:
					case HCLOCK_DUMMY+112:
					case HCLOCK_DUMMY+114:
					case HCLOCK_DUMMY+116:
					case HCLOCK_DUMMY+118:
					case HCLOCK_DUMMY+120:
					case HCLOCK_DUMMY+122:
					case HCLOCK_DUMMY+124:
					case HCLOCK_DUMMY+126:
					case HCLOCK_DUMMY+128:
					case HCLOCK_DUMMY+130:
					case HCLOCK_DUMMY+132:
					case HCLOCK_DUMMY+134:
					case HCLOCK_DUMMY+136:
					case HCLOCK_DUMMY+138:
					case HCLOCK_DUMMY+140:
					case HCLOCK_DUMMY+142:
					case HCLOCK_DUMMY+144:
					case HCLOCK_DUMMY+146:
					case HCLOCK_DUMMY+148:
					case HCLOCK_DUMMY+150:
					case HCLOCK_DUMMY+152:
					case HCLOCK_DUMMY+154:
					case HCLOCK_DUMMY+156:
					case HCLOCK_DUMMY+158:
					case HCLOCK_DUMMY+160:
					case HCLOCK_DUMMY+162:
					case HCLOCK_DUMMY+164:
					case HCLOCK_DUMMY+166:
					case HCLOCK_DUMMY+168:
					case HCLOCK_DUMMY+170:
					case HCLOCK_DUMMY+172:
					case HCLOCK_DUMMY+174:
					case HCLOCK_DUMMY+176:
					case HCLOCK_DUMMY+178:
					case HCLOCK_DUMMY+180:
					case HCLOCK_DUMMY+182:
					case HCLOCK_DUMMY+184:
					case HCLOCK_DUMMY+186:
					case HCLOCK_DUMMY+188:
					case HCLOCK_DUMMY+190:
					case HCLOCK_DUMMY+192:
					case HCLOCK_DUMMY+194:
					case HCLOCK_DUMMY+196:
					case HCLOCK_DUMMY+198:
					case HCLOCK_DUMMY+200:
					case HCLOCK_DUMMY+202:
					case HCLOCK_DUMMY+204:
					case HCLOCK_DUMMY+206:
					case HCLOCK_DUMMY+208:
					case HCLOCK_DUMMY+210:
					case HCLOCK_DUMMY+212:
					case HCLOCK_DUMMY+214:
					case HCLOCK_DUMMY+216:
					case HCLOCK_DUMMY+218:
					case HCLOCK_DUMMY+220:
					case HCLOCK_DUMMY+222:
					case HCLOCK_DUMMY+224:
					case HCLOCK_DUMMY+226:
					case HCLOCK_DUMMY+228:
					case HCLOCK_DUMMY+230:
					case HCLOCK_DUMMY+232:
					case HCLOCK_DUMMY+234:
					case HCLOCK_DUMMY+236:
					case HCLOCK_DUMMY+238:
					case HCLOCK_DUMMY+240:
					case HCLOCK_DUMMY+242:
					case HCLOCK_DUMMY+244:
					case HCLOCK_DUMMY+246:
					case HCLOCK_DUMMY+248:
					case HCLOCK_DUMMY+250:
					case HCLOCK_DUMMY+252:
					case HCLOCK_DUMMY+254:
					case HCLOCK_DUMMY+256:
					case HCLOCK_DUMMY+258:
					case HCLOCK_DUMMY+260:
					case HCLOCK_DUMMY+262:
					case HCLOCK_DUMMY+264:
					case HCLOCK_DUMMY+266:
					case HCLOCK_DUMMY+268:
					case HCLOCK_DUMMY+270:
					case HCLOCK_DUMMY+272:
					case HCLOCK_DUMMY+274:
					case HCLOCK_DUMMY+276:
					case HCLOCK_DUMMY+278:
					case HCLOCK_DUMMY+280:
					case HCLOCK_DUMMY+282:
					case HCLOCK_DUMMY+284:
					case HCLOCK_DUMMY+286:
					case HCLOCK_DUMMY+288:
					case HCLOCK_DUMMY+290:
					case HCLOCK_DUMMY+292:
					case HCLOCK_DUMMY+294:
					case HCLOCK_DUMMY+296:
					case HCLOCK_DUMMY+298:
					case HCLOCK_DUMMY+300:
					case HCLOCK_DUMMY+302:
					case HCLOCK_DUMMY+304:
					case HCLOCK_DUMMY+306:
					case HCLOCK_DUMMY+308:
					case HCLOCK_DUMMY+310:
					case HCLOCK_DUMMY+312:
					case HCLOCK_DUMMY+314:
					case HCLOCK_DUMMY+316:
					{
						NST_COMPILE_ASSERT( HCLOCK_DUMMY & 1 );

						cycles.hClock = cycles.count | 1;

						if (cycles.count <= HCLOCK_DUMMY+318)
							break;
					}

					case HCLOCK_DUMMY+318:

						cycles.hClock = 320;
						cycles.vClock += HCLOCK_DUMMY;
						cycles.count -= HCLOCK_DUMMY;

						if (cycles.count <= 320)
							break;

						goto HBlankOff;

					default:

						NST_UNREACHABLE();
				}
			}

			cycles.count = GetCycles();
		}
	}
}
