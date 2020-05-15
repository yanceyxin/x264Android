#!/bin/sh

ANDROID_NDK=$HOME/Library/Android/sdk//android-ndk-r14b
SYSROOT=$ANDROID_NDK/platforms/android-19/arch-arm
CROSS_PREFIX=$ANDROID_NDK/toolchains/arm-linux-androideabi-4.9/prebuilt/darwin-x86_64/bin/arm-linux-androideabi-
EXTRA_CFLAGS="-march=armv7-a -mfloat-abi=softfp -mfpu=neon -D__ANDROID__ -D__ARM_ARCH_7__ -D__ARM_ARCH_7A__"
EXTRA_LDFLAGS="-nostdlib"
PREFIX=`pwd`/libs/armeabi-v7a
echo $CROSS_PREFIX
echo $SYSROOT
echo $EXTRA_CFLAGS
echo $EXTRA_LDFLAGS

./configure --prefix=$PREFIX \
        --host=arm-linux \
        --sysroot=$SYSROOT \
        --cross-prefix=$CROSS_PREFIX \
        --extra-cflags="$EXTRA_CFLAGS" \
        --extra-ldflags="$EXTRA_LDFLAGS" \
        --bit-depth=8 \
        --chroma-format=420 \
        --disable-interlaced \
        --enable-pic \
        --enable-static \
        --enable-strip \
        --disable-cli \
        --disable-win32thread \
        --disable-avs \
        --disable-swscale \
        --disable-lavf \
        --disable-ffms \
        --disable-gpac \
        --disable-lsmash \
        --disable-opencl

make clean
make STRIP= -j8 install || exit 1
