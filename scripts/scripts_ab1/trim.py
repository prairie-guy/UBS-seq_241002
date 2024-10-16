#!/usr/bin/env python3

def trim(seq, cutoff=0.05):
        """Trims the sequence using Richard Mott's modified trimming algorithm.

        Keyword argument:
        seq -- sequence to be trimmed
        cutoff -- probability cutoff value

        Trimmed bases are determined from their segment score, ultimately
        determined from each base's quality values.

        More on:
        http://www.phrap.org/phredphrap/phred.html
        http://www.clcbio.com/manual/genomics/Quality_trimming.html
        """
        # set flag for trimming
        start = False
        # set minimum segment size
        segment = 20
        trim_start = 0

        if len(seq) <= segment:
            raise ValueError('Sequence can not be trimmed because \
                             it is shorter than the trim segment size')
        else:
            # calculate probability back from formula used
            # to calculate phred qual values
            score_list = [cutoff - (10 ** (qual/-10.0)) for
                         qual in self.qual_val]

            # calculate cummulative score_list
            # if cummulative value < 0, set to 0
            # first value is set to 0 (assumption: trim_start is always > 0)
            running_sum = [0]
            for i in range(1, len(score_list)):
                num = running_sum[-1] + score_list[i]
                if num < 0:
                    running_sum.append(0)
                else:
                    running_sum.append(num)
                    if not start:
                        # trim_start = value when cummulative starts to be > 0
                        trim_start = i
                        start = True

            # trim_finish = index of the highest cummulative value,
            # marking the segment with the highest cummulative score
            trim_finish = running_sum.index(max(running_sum))

            return seq[trim_start:trim_finish]
