#!/usr/bin/env nu

# Key = `SeqID`
# Join an experiment.csv with the correposnding sequence data with the key `SeqID`
# `experiment.csv` is assumed to have a column `name` with each containing `BAT-nn` and `R1|R2`
# Experimental Protocol should be separated by two sets of -> -> Within, each experimental condition should be comma-delimited

def get-sequences  [dir] {
    let fns = ls $dir
    ls $dir
    | get name
    #| parse --regex '.*(?P<bat>BAT-\d+).*(?P<std>R[1|2]).*'
    #| parse --regex '.*(?P<bat>BAT-\d+).*(?P<std>R[1|2]).*'
    | merge $fns 
    | rename SeqID Std Sequence
    | select SeqID Std Sequence
    }

def create-samples [exp_path_csv seq_dir] {
    let seqs = get-sequences $seq_dir
    open $exp_path_csv
    | select SampleID SeqID Content
    | join $seqs SeqID
    | insert FullID {|row| $'($row.SampleID)_($row.Std)'}
    | insert Experiment {|row| $row.Content | str replace -r '.*->\W*(.*)\W*->.*' '$1'}
    | update Experiment {|row| $row.Experiment | str replace --all 'for' ';'}
    | update Experiment {|row| $row.Experiment | str replace --all ';' '; '}
    | select SeqID SampleID FullID Experiment Std Sequence
    | sort-by SeqID -n
    }



def main [exp_path_csv seq_dir] {
    create-samples $exp_path_csv $seq_dir
    |to csv
    |save -f 'samples.csv'}

# |parse --regex '.*->\W(?P<desc>.*)\W->.*'
