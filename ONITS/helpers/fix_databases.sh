# This script fixes the trimmed databases

# TODO:check entries
euk=$1
unite=$2

sed -i'.bak' 's/tax=//' "$euk"
sed -i'.bak' 's/,$//' "$euk"

sed -i'.bak' 's/;tax=/|/' "$unite"
sed -i'.bak' 's/,tax=/;tax=/' "$unite"
sed -i'.bak' 's/,$//' "$unite"
