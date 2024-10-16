#!/usr/bin/env nu

# Compare two converted csv files, returning  average converted 'c' and 'm5c' for each
def compare-converted [converted_1_csv: path converted_2_csv: path]: nothing -> table  {
    let c1 = open $converted_1_csv | rename ref sample ac_1 a5c_1
    let c2 = open $converted_2_csv | select ac a5c | rename ac_2 a5c_2
    $c1 | merge $c2 | select ref sample ac_1 ac_2 a5c_1 a5c_2
}

def main [c1 c2] {
    compare-converted $c1 $c2 | save -f avg_converted_chang_and_bryan.csv
}
