BEGIN{FS="\t"; OFS="/t"}
{print $1"_"$2"_"$3"_"$4"_.""\t"$5}
