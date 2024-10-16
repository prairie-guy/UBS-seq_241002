#!/usr/bin/env nu

# Key = `SampeID`
# Join an experiment.csv with the correposnding sequence data with the key `SampleID`
# `experiment.csv` is assumed to have a column `name` with each containing `BAT-nn` and `R1|R2`
# Experimental Protocol should be separated by two sets of -> -> Within, each experimental condition should be comma-delimited

def get-sequences [dir] {
    let fns = ls $dir | where name =~ 'gz$' and name !~ 'tar\.gz$'
    $fns
    |get name
    | parse --regex '.*-(?P<SeqID>.+)_.*_(?P<Std>R[1|2]).*'
    |merge $fns
    |rename SampleID Std Sequence
    |select SampleID Std Sequence
    |update SampleID  {|it|if $it.SampleID == '0' {'O'} else {$it.SampleID}}
    |sort-by SampleID -n
    }

def create-samples [exp_path_csv seq_dir] {
    let seqs = get-sequences $seq_dir
    open $exp_path_csv
    | select SampleID  Content
    | join $seqs SampleID
    | insert FullID {|row| $'($row.SampleID)_($row.Std)'}
    | insert Experiment {|row| $row.Content | str replace -r '.*->\W*(.*)\W*->.*' '$1'}
    | update Experiment {|row| $row.Experiment | str replace --all 'for' ';'}
    | update Experiment {|row| $row.Experiment | str replace --all ';' '; '}
    | select SampleID FullID Experiment Std Sequence
    | sort-by SampleID -n
    }


def main [exp_path_csv seq_dir] {
    create-samples $exp_path_csv $seq_dir
    |to csv
    |save -f 'samples.csv'}

#|parse -r '.*/(?P<SampleID>\w\d*)_.*(?P<Std>R[12]).*'
