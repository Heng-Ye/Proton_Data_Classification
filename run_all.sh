#!/bin/bash

class_name="ProtonNewTreeMaker_run"
class_namex=$class_name"X"
echo $class_namex


RUN[0]=5219
RUN[1]=5225
RUN[2]=5235
RUN[3]=5240
RUN[4]=5244
RUN[5]=5308
RUN[6]=5311
RUN[7]=5315
RUN[8]=5338
RUN[9]=5387
RUN[10]=5423
RUN[11]=5424
RUN[12]=5426
RUN[13]=5455
RUN[14]=5456
RUN[15]=5457
RUN[16]=5458
RUN[17]=5460
RUN[18]=5809
RUN[19]=5810
RUN[20]=5814
RUN[21]=5816
RUN[22]=5817
RUN[23]=5842
RUN[24]=5843
RUN[25]=5844

for run in ${RUN[*]}
do
    printf "%s\n" $run

    selector_name=$class_name$run
    ana_name="makeproton_ana_"$class_name$run

    header_name=$class_name$run".h"
    tmp_header_name=$class_name${run}"_tmp.h"

    class_code=$class_name${run}".C"
    tmp_class_code=$class_name$run"_tmp.C"

    g++ makeproton_ana.cc `root-config --libs --cflags` -o $ana_name
    ./$ana_name "./file_list/file_run"${run}"_reco2.txt" $selector_name
    sed '/Init(tree)\;/i if (tree-\>InheritsFrom(\"TChain\")) ((TChain\*)tree)-\>LoadTree(0);' $header_name > $tmp_header_name
    mv $tmp_header_name $header_name
    cp -prv $class_namex".C"  $tmp_class_code
    sed 's/5387/'${run}'/g' $tmp_class_code > $class_code
    rm -f $tmp_class_code
    root_exe_str="root -b -q 'RunAna.C(\""$selector_name\"")'"

    echo $root_exe_str" ......"
    eval $root_exe_str

done


