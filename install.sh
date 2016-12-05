# Download latest version
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip

#install the database
java -jar snpEff/snpEff.jar download GRCh38.86
