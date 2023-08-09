for i in *_platon.csv ; do
    platon=$i
    blast=${i%"_platon.csv"}"_blast.csv"

    i_platon=`cat $platon | cut -f6 -d "," | tail -n +2`
    i_platon="${i_platon//$'\n'/ }"
    i_blast=`cat $blast | cut -f2 -d "," | tail -n +2`
    i_blast="${i_blast//$'\n'/ }"
    for item1 in ${i_platon[@]}; do
        for item2 in ${i_blast[@]}; do
            if [[ "$item1" == "$item2" ]]; then
                intersections+=( "$item1" )
                break
            fi
        done
    done
    echo -e "${intersections[@]}" > ${i%"_platon.csv"}"_common.txt"
    unset intersections[@]
done
