# Work log

## Successful commands
- `date +%y-%m-%d_%H-%M && echo $$`
- `CPPFLAGS='-I/opt/sw/include' CFLAGS='-I/opt/sw/include' LDFLAGS='-L/opt/sw/lib' ./configure --disable-libcurl --with-libdeflate --enable-plugins`
- `perl -0pi -e 's/-ldeflate/-l:libdeflate.a/g; s/-lbz2/-l:libbz2.a/g; s/-llzma/-l:liblzma.a/g' config.mk && rg -n '^LIBS|^LDFLAGS' config.mk`
- `make -j4 lib-static bgzip tabix`
- `rsync -av htslib /opt/sw/include/ && rsync -av libhts.a /opt/sw/lib/ && rsync -av bgzip tabix /opt/sw/bin/`
- `wget -O ../bcftools-1.23.1.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2 && tar -xjf ../bcftools-1.23.1.tar.bz2 -C ..`
- `for f in score.c score.h munge.c liftover.c metal.c blup.c pgs.c pgs.mk; do wget -O "plugins/$f" "https://raw.githubusercontent.com/freeseek/score/master/$f"; done`
- `wget -O assoc_plot.R https://raw.githubusercontent.com/freeseek/score/master/assoc_plot.R`
- `CPPFLAGS='-I/opt/sw/include -I/usr/include/suitesparse' CFLAGS='-I/opt/sw/include -I/usr/include/suitesparse' LDFLAGS='-L/opt/sw/lib -L/usr/lib/x86_64-linux-gnu' ./configure --enable-plugins --prefix=/opt/sw --with-htslib=../htslib-1.23.1 --with-bcf-plugin-dir=/opt/sw/bcf-plugins`
- `make htslib_static.mk`
- `perl -0pi -e 's/^(HTSLIB_LIB\s*=).*/$1 -l:libhts.a -l:libdeflate.a -l:libbz2.a -l:liblzma.a -lz -lm -ldl -lpthread/m; s/^(HTSLIB_LDFLAGS\s*=).*/$1 -L..\/htslib-1.23.1 -L\/opt\/sw\/lib/m' config.mk`
- `make -j4`
- `make install && install -m 0755 assoc_plot.R /opt/sw/bin/assoc_plot.R`
- `/opt/sw/bin/bcftools --version`
- `BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins /opt/sw/bin/bcftools plugin -lv | rg '^-- (score|munge|liftover|metal|blup|pgs) --'`
- `ldd /opt/sw/bin/bcftools | rg 'libhts|libdeflate|libbz2|liblzma' || true`
- `ldd /opt/sw/bcf-plugins/pgs.so | rg 'cholmod|amd|colamd|camd|ccolamd|suitesparse|blas|lapack|gomp|gfortran'`
- `for p in score munge liftover metal blup pgs; do echo "== $p =="; BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins /opt/sw/bin/bcftools +$p -h >/tmp/bcfplugin-$p.help 2>&1; echo rc=$?; head -8 /tmp/bcfplugin-$p.help; done`
- `file /opt/sw/bin/bcftools /opt/sw/bcf-plugins/pgs.so /opt/sw/bin/assoc_plot.R && ls -l /opt/sw/bin/assoc_plot.R /opt/sw/bcf-plugins/{score,munge,liftover,metal,blup,pgs}.so`
- `ldd /opt/sw/bin/bcftools`

## Successful actions
- Wrote the approved plan to `plan_26-04-28_18-07_1808206.md`.
- Built htslib 1.23.1 with static libdeflate, bzip2, and lzma linkage.
- Installed htslib headers, `libhts.a`, `bgzip`, and `tabix` into `/opt/sw`.
- Downloaded and extracted bcftools 1.23.1.
- Added current freescore plugin sources from `freeseek/score`.
- Configured bcftools to use `../htslib-1.23.1` and `/opt/sw/bcf-plugins`.
- Patched bcftools linkage to use static `libhts.a`, `libdeflate.a`, `libbz2.a`, and `liblzma.a`.
- Attempted static SuiteSparse/CHOLMOD linkage for `pgs.so`; changed only `pgs.so` to dynamic system CHOLMOD after non-PIC SuiteSparse/METIS archive objects prevented static plugin linkage.
- Built and installed bcftools 1.23.1 and plugins into `/opt/sw`.
- Installed `assoc_plot.R` into `/opt/sw/bin`.
- Verified bcftools 1.23.1 and htslib 1.23.1 runtime versions.
- Verified freescore plugin discovery for `score`, `munge`, `liftover`, `metal`, `blup`, and `pgs`.
- Verified the main bcftools binary has no dynamic `libhts`, `libdeflate`, `libbz2`, or `liblzma` dependencies.
- Verified `pgs.so` uses dynamic system CHOLMOD/SuiteSparse/OpenBLAS dependencies.
