#!/bin/bash

fastq=$1
prefix=$2
startk=$3


if [ $# -lt 3 ]
then
	echo -e "The script takes three mandatory arguments: 1) input fastq file 2) prefix for output files 3) start k-mer length"
	echo -e "Example: ./bless_iterate_over_ks.sh seqs.fastq.gz prefix 21"
	echo -e "\nOptional one can also specify the maximum memory usage (parameter 4 [GB]) and the number of threads to be used (parameter 5 [int])"
	echo -e "a 6th argument ('clean' in the example below - could be anything though) will trigger removal of the output fastq files,"
	echo -e "i.e. sometimes one will only be interested in the number errors that were corrected"
	echo -e "Example (keep fastq files): ./bless_iterate_over_ks.sh seqs.fastq.gz prefix 21 10 10"
	echo -e "Example (remove fastq files): ./bless_iterate_over_ks.sh seqs.fastq.gz prefix 21 10 10 clean"
	exit 1
fi


#######

endk=$(( startk + 60 ))
kmercountok=0
ok=0
k=$startk
while [ $ok -eq 0 ]
do
	if [ ! -f "bless-k$k.log" ]
	then
		echo -e "\n$(date)\tRunning bless with k=$k\n"
#		/src/bless/v1p02/bless -read $fastq -kmerlength $k -notrim -prefix $prefix-k$k > bless-k$k.log
		cmd="bless -read $fastq -kmerlength $k -notrim -prefix $prefix-k$k"
		if [ $# -gt 3 ]
		then
			cmd="$cmd -max_mem $4 -smpthread $5"
		fi
		if [ ! $6 ]
		then
			cmd="$cmd -gzip"
		fi
		echo -e "$cmd"
		#in clean mode write corrected reads to /dev/null, so not to use up space
		if [ $6 ]
		then
			ln -s /dev/null $prefix-k$k.corrected.fastq.00000
		fi
		$cmd > bless-k$k.log.tmp
		#when finished, rename the log so that it's obvious if a process stopped before completed. Only complete ones will be evaluated subsequently
		mv bless-k$k.log.tmp bless-k$k.log
		#in clean mode remove empty read file
		if [ $6 ]
		then
			rm $prefix-k$k.corrected.fastq
		fi
	fi

	if [ "$kmercountok" -eq 0 ]
	then
		kmercount=$(cat bless-k$k.log | grep "Number of unique solid k-mers" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"')
		echo -e "$(date)\tNumber of unique solid k-mers (k=$k): Ns = $kmercount"
		#check if number of kmers/4^k < 0.0001
		out=$(echo $(cat bless-k$k.log | grep "k-mer length" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"') $kmercount | perl -ne 'chomp; @a=split(" "); $res=($a[1]/(4**$a[0])); if ($res < 0.0001){print "$res -> OK!\n"}else{print "$res -> NOT OK!\n"}')
		echo -e "$(date)\tChecking criterion: Ns / 4 ^ k < 0.0001 .."
		echo -e "$(date)\t$kmercount / 4 ^ $k = $out"
		if [ "$(echo "$out" | grep "NOT" |wc -l)" -eq "0" ]
		then
			kmercountok=1
			echo -e "$(date)\tWill now search for the k-mer length at which the number of corrected bases is maximal"
		else
			if [ "$(find ./ -name "$prefix-k$k*")" ]
			then
				echo -e "$(date)\tremoving results for k=$k"
				rm -v $prefix-k$k*
			fi
		fi
		echo -e "$(date)\tmoving on to next k"
		
	else
		first="bless-k$(( k - 2 )).log"
		second="bless-k$k.log"
#		first=$(ls -hrlt bless-k* | tail -n 2 | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"'| head -n 1)
#		second=$(ls -hrlt bless-k* | tail -n 2 | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"'| tail -n 1)
		#echo $first $second
		firstk=$(cat $first | grep "k-mer length" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"')
		secondk=$(cat $second | grep "k-mer length" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"')
		#echo $firstk $secondk
		firstcount=$(cat $first | grep "Number of corrected errors" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"')
		secondcount=$(cat $second | grep "Number of corrected errors" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"')
		#echo $firstcount $secondcount

		echo -e "$(date)\tChecking number of corrected errors - k=$firstk vs. k=$secondk"
		
		if [ "$firstcount" -gt "$secondcount" ]
		then
			echo -e "$(date)\t$firstk is the winner - $firstcount > $secondcount"
			echo -e "$(date)\tDone!\n"
			rm $prefix-k$secondk*
			echo -e "$(( k - 2 ))" > $prefix.bestk
			ok=1
		else
			echo -e "$(date)\t$secondk is the winner - $firstcount < $secondcount"
			echo -e "$(date)\tmoving on to next k"
			if [ "$(find ./ -name "$prefix-k$firstk*")" ]
			then
				echo -e "$(date)\tremoving results for k=$firstk"
				rm $prefix-k$firstk*
			fi
		fi

	fi

	if [ "$(( k + 2 ))" -gt "$endk" ]
	then
		ok=1
		echo -e "Reached upper k-mer limit ($endk)\n"
	fi

	k=$(( k + 2 ))

done

