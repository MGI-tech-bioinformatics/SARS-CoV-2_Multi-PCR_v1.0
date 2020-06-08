#!/bin/bash

workdir=`pwd`

mkdir -p $workdir/tools/packages

SEQTK_URL="https://github.com/lh3/seqtk.git"
SOAPnuke_URL="https://github.com/BGI-flexlab/SOAPnuke/archive/1.5.6-linux.tar.gz"
bamdst_URL="https://github.com/shiquan/bamdst.git"
bcftools_URL="https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2"
tabix_URL="https://github.com/samtools/tabix.git"

RM_TMP()
{
    local FILE=$1

    if [ -f $FILE ] || [ -d $FILE ]
    then
        rm -rf $FILE
    fi
}

echo "Installing seqtk ..."
cd $workdir/tools/packages
RM_TMP $workdir/tools/seqtk
RM_TMP $workdir/tools/packages/seqtk
git clone $SEQTK_URL
cd seqtk; make
cd $workdir/tools; ln -s ./packages/seqtk/seqtk
echo "Installing seqtk successfully"

echo "Installing SOAPnuke ..."
cd $workdir/tools/packages
RM_TMP $workdir/tools/SOAPnuke
RM_TMP $workdir/tools/packages/SOAPnuke-1.5.6-linux
RM_TMP $workdir/tools/packages/1.5.6-linux.tar.gz
wget -c $SOAPnuke_URL
tar -zxvf 1.5.6-linux.tar.gz
cd $workdir/tools; ln -s ./packages/SOAPnuke-1.5.6-linux/bin/SOAPnuke
echo "Installing SOAPnuke successfully"

echo "Installing bamdst ..."
cd $workdir/tools/packages
RM_TMP $workdir/tools/bamdst
RM_TMP $workdir/tools/packages/bamdst
git clone $bamdst_URL
cd bamdst; make
cd $workdir/tools; ln -s ./packages/bamdst/bamdst
echo "Installing bamdst successfully"

echo "Installing bcftools ..."
cd $workdir/tools/packages
RM_TMP $workdir/tools/bcftools
RM_TMP $workdir/tools/packages/bcftools
RM_TMP $workdir/tools/packages/bcftools-1.8.tar.bz2
wget -c $bcftools_URL
tar jxf bcftools-1.8.tar.bz2
cd bcftools-1.8; ./configure; make
cd $workdir/tools; ln -s ./packages/bcftools-1.8/bcftools
echo "Installing bcftools successfully"

echo "Installing tabix and bgzip ..."
cd $workdir/tools/packages
RM_TMP $workdir/tools/tabix
RM_TMP $workdir/tools/bgzip
RM_TMP $workdir/tools/packages/tabix
git clone $tabix_URL
cd tabix; make
cd $workdir/tools; ln -s ./packages/tabix/tabix; ln -s ./packages/tabix/bgzip
echo "Installing tabix and bgzip successfully"
