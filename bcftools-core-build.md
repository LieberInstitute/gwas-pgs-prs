Goal: 
 -build and install bcftools with the plugins from https://github.com/freeseek/score
 -make sure libdeflate is used, as installed under $SW prefix, 
  source tree and build available in ../libdeflate
 -use htslib already downloaded in ../htslib-* (highest version available, currently 1.23.1)


## Software installation prefix

Use SW=/opt/sw as the installation prefix;

Download the htslib-matching version of bcftools from https://www.htslib.org/download/ 
and build it in this directory, with the plugins and all the auxiliary files as described at

https://github.com/freeseek/score#installation
and
https://software.broadinstitute.org/software/score/

## Semi-static build

I want to make sure that the bcftools version is linked statically with htslib.a and libdeflate.a 
(not dynamically with the .so files).

The user may ask to also link statically bzip2 and lzma libraries. 

The recipe the user used in the past for that, assumed we have libdeflate, bzip2 and lzma static 
libraries installed under $SW prefix. 
(If not, look for their source/build folders in the parent directory, or build them there)
---
LDFLAGS="-L$SW/lib"
CPPFLAGS="-I$SW/include"
CFLAGS="$CPPFLAGS"
./configure --disable-libcurl --with-libdeflate --enable-plugins
perl -i -pe 's/-ldeflate/-l:libdeflate.a/;s/-lbz2/-l:libbz2.a/;s/-llzma/-l:liblzma.a/' config.mk
make -j4 lib-static
make bgzip tabix
## optional, to make only the static library available:
rsync -av htslib $SW/include/ 
rsync -av libhts.a $SW/lib/
rsync -av bgzip tabix $SW/bin/
---

Then bctools can be built similarly and linked against this htslib version, like this:

export BCFTOOLS_PLUGINS=$SW/bcf-plugins
./configure --enable-plugins --prefix=$SW --with-htslib=../htslib-1.23.1 --with-bcf-plugin-dir=$BCFTOOLS_PLUGINS
perl -i -pe 's/^(HTSLIB_LIB\s*=).+/$1 -l:libhts.a -lcurl -lcrypto -l:libdeflate.a -l:libbz2.a -l:liblzma.a -lm -ldl/' config.mk
make -j4 && make install
