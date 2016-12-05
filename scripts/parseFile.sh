VCF=$1
cat $VCF | awk 'grep '0|0' $5'
