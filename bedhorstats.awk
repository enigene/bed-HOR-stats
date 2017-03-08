# Collect stats from BED files after file(s) has been passed through hmmertblout2bed.awk
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
# v1.1, 19 Dec 2016 - Add filename stamp
# v1.0, 10 Nov 2016 - Initial release
# 
# Usage: awk -v l="AS-SF-HORs-list-s.txt" -v scorestat=1 -v lengthstat=1 -f bedhorstats-v1.awk input.bed > output.stats

BEGIN {
# HOR list preparation
    if (hNum == 0) {
    	while ((getline __input < l) > 0) {
    # Add all HORs from list
        	if (__input ~ /^S/) {
        		fullName = __input;
    			dotpos = index(fullName, ".");
    			nameBeforeDot = substr(fullName, 1, dotpos-1);

    			if (!(nameBeforeDot in horList)) {
    				horList[++hNum] = nameBeforeDot;
    				horList[nameBeforeDot] = hNum;
    			}

        		horMonList[++mNum] = fullName;
        		horMonList[fullName] = mNum;
        	}
    # Add HOR groups
            if (__input ~ /^g/) {
                fullName = __input;
                split(fullName, tmpA);
                gNum = substr(tmpA[1], 2);
                gNums[gNum]++;
                gHorName = tmpA[2];
                horComment = tmpA[3];
                horGroups[gNum][++gHorNum] = gHorName;
                horGroups[gNum][gHorName] = gHorNum;
                horComments[gNum][gHorName] = horComment;
            }
        }
        close(l);
    }

    first = 1;
}

NR == 1 {
  fbase = FILENAME;
  fext = FILENAME;
  gsub(/.*[/]/, "", fbase);
  sub(/[.][^.]+$/, "", fbase);
  sub(/.*[.]/, "", fext);
}

{
	fullName = $4;
	dotpos = index(fullName, ".");
	nameBeforeDot = substr(fullName, 1, dotpos-1);
	start = $2;
	finish = $3;
	score = sprintf("%.2f",$5);
	monLength = finish - start;

	horA[nameBeforeDot]++;
	monA[fullName]++;
	horScoreA[nameBeforeDot][score]++;
	horLengthA[nameBeforeDot][monLength]++;
	monScoreA[fullName][score]++;
	monLengthA[fullName][monLength]++;

	scA[sprintf("%d",score)];

	if (first) {
		maxScore = minScore = score;
		maxLength = minLength = monLength;
		first = 0;
	}
}

END {
# Print elements in SF
	print "#" fbase;
	print "#HOR_list";
	if (scorestat != 0) {
		printf("#%s\t\t\n", fbase) > fbase "-score-stat.out";
		printf("%s\n", "#Distribution_of_scores") >> fbase "-score-stat.out";
		printf("%s\t%s\t%s\t%s\n", "HOR_list", "S<=61", "61<S<97", "S>=97") >> fbase "-score-stat.out";
	}
	if (lengthstat != 0) {
		printf("#%s\t\t\n", fbase) > fbase "-length-stat.out";
		printf("%s\n", "#Distribution_of_lengths") >> fbase "-length-stat.out";
		printf("%s\t%s\t%s\t%s\n", "HOR_list", "L<=50", "50<L<76", "L>=76") >> fbase "-length-stat.out";
	}
	for (iInput in horA) {
		if (iInput ~ /^S1/) sf1c += horA[iInput];
		if (iInput ~ /^S2/) sf2c += horA[iInput];
		if (iInput ~ /^S3/) sf3c += horA[iInput];
		if (iInput ~ /^S5/) sf5c += horA[iInput];
		if (iInput ~ /^S6/) sf6c += horA[iInput];
#		if (iInput ~ /^S4\/6/) sf4_6c += horA[iInput];
		if (iInput ~ /^S4/) sf4c += horA[iInput];
	}
	printf("%s\t%d\n", "SF1", sf1c+0);
	printf("%s\t%d\n", "SF2", sf2c+0);
	printf("%s\t%d\n", "SF3", sf3c+0);
	printf("%s\t%d\n", "SF5", sf5c+0);
	printf("%s\t%d\n", "SF6", sf6c+0);
	printf("%s\t%d\n", "SF4", sf4c+0);

# Print elements in HOR list
	for (iList=1; iList<=hNum; iList++) {
		if (horList[iList] in horA) {
			if (horA[horList[iList]]+0 > 0) {
			  hCount++;
			}
			printf("%s\t%d\n", horList[iList], horA[horList[iList]]+0);
# Score distribution stats
			if (scorestat != 0) {
				for (sval in horScoreA[horList[iList]]) {
					if (sval <= 61) {
						horScoreDistrA[horList[iList]]["leftCount"] += horScoreA[horList[iList]][sval];
						totalHorScoreA["leftCount"] += horScoreA[horList[iList]][sval];
					}
					if ((sval > 61)&&(sval < 97)) {
						horScoreDistrA[horList[iList]]["middleCount"] += horScoreA[horList[iList]][sval];
						totalHorScoreA["middleCount"] += horScoreA[horList[iList]][sval];
					}
					if (sval >= 97) {
						horScoreDistrA[horList[iList]]["rightCount"] += horScoreA[horList[iList]][sval];
						totalHorScoreA["rightCount"] += horScoreA[horList[iList]][sval];
					}
#					printf("%s\t%.2f\t%d\n", horList[iList], sval, horScoreA[horList[iList]][sval]) >> fbase "-score-stat.out";
				}
				printf("%s\t%d\t%d\t%d\n", horList[iList], \
					horScoreDistrA[horList[iList]]["leftCount"], \
					horScoreDistrA[horList[iList]]["middleCount"], \
					horScoreDistrA[horList[iList]]["rightCount"]) >> fbase "-score-stat.out";
			}
# Length distribution stats
			if (lengthstat != 0) {
				for (lval in horLengthA[horList[iList]]) {
					if (lval <= 50) {
						horLengthDistrA[horList[iList]]["leftCount"] += horLengthA[horList[iList]][lval];
						totalHorLengthA["leftCount"] += horLengthA[horList[iList]][lval];
					}
					if ((lval > 50)&&(lval < 76)) {
						horLengthDistrA[horList[iList]]["middleCount"] += horLengthA[horList[iList]][lval];
						totalHorLengthA["middleCount"] += horLengthA[horList[iList]][lval];
					}
					if (lval >= 76) {
						horLengthDistrA[horList[iList]]["rightCount"] += horLengthA[horList[iList]][lval];
						totalHorLengthA["rightCount"] += horLengthA[horList[iList]][lval];
					}
#					printf("%s\t%.2f\t%d\n", horList[iList], lval, horLengthA[horList[iList]][lval]) >> fbase "-length-stat.out";
				}
				printf("%s\t%d\t%d\t%d\n", horList[iList], \
					horLengthDistrA[horList[iList]]["leftCount"], \
					horLengthDistrA[horList[iList]]["middleCount"], \
					horLengthDistrA[horList[iList]]["rightCount"]) >> fbase "-length-stat.out";
			}
		} else {
			printf("%s\t%d\n", horList[iList], horA[horList[iList]]+0);
			if (scorestat != 0) {
				printf("%s\t%d\t%d\t%d\n", horList[iList], \
					horScoreDistrA[horList[iList]]["leftCount"], \
					horScoreDistrA[horList[iList]]["middleCount"], \
					horScoreDistrA[horList[iList]]["rightCount"]) >> fbase "-score-stat.out";
			}
			if (lengthstat != 0) {
				printf("%s\t%d\t%d\t%d\n", horList[iList], \
					horLengthDistrA[horList[iList]]["leftCount"], \
					horLengthDistrA[horList[iList]]["middleCount"], \
					horLengthDistrA[horList[iList]]["rightCount"]) >> fbase "-length-stat.out";
			}
		}
	}

# Print total for scores
	if (scorestat != 0) {
		printf("%s\t%d\t%d\t%d\n", "Total_in_group:", \
			totalHorScoreA["leftCount"], \
			totalHorScoreA["middleCount"], \
			totalHorScoreA["rightCount"]) >> fbase "-score-stat.out";
	}
# Print total for lengts
	if (lengthstat != 0) {
		printf("%s\t%d\t%d\t%d\n", "Total_in_groups:", \
			totalHorLengthA["leftCount"], \
			totalHorLengthA["middleCount"], \
			totalHorLengthA["rightCount"]) >> fbase "-length-stat.out";
	}

# Print elements absent in HOR list
	for (iInput in horA) {
		if (!(iInput in horList)) {
			foundAbsentHORs++;
			if (foundAbsentHORs == 1) {
				print "\n#Elements_absent_in_HOR_list:";
			}
			printf("%s\t%d\n", iInput, horA[iInput]+0);
		}
	}

	printf("%s\t%d\n\n", "HORs_with_non-zero_amount_of_monomers:", hCount+0);
	
	print "#List_of_HOR_monomers";
	if (scorestat != 0) {
		printf("\n%s\t%s\t%s\t%s\n", "List_of_HOR_monomers", "S<=61", "61<S<97", "S>=97") >> fbase "-score-stat.out";
	}
	if (lengthstat != 0) {
		printf("\n%s\t%s\t%s\t%s\n", "List_of_HOR_monomers", "L<=50", "50<L<76", "L>=76") >> fbase "-length-stat.out";
	}

# Print elements in HOR mon list
	for (iList=1; iList<=mNum; iList++) {
		if (horMonList[iList] in monA) {
			mCount += monA[horMonList[iList]];
			printf("%s\t%d\n", horMonList[iList], monA[horMonList[iList]]+0);
# Score distribution stats
			if (scorestat != 0) {
				for (sval in monScoreA[horMonList[iList]]) {
					if (sval <= 61) {
						monScoreDistrA[horMonList[iList]]["leftCount"] += monScoreA[horMonList[iList]][sval];
					}
					if ((sval > 61)&&(sval < 97)) {
						monScoreDistrA[horMonList[iList]]["middleCount"] += monScoreA[horMonList[iList]][sval];
					}
					if (sval >= 97) {
						monScoreDistrA[horMonList[iList]]["rightCount"] += monScoreA[horMonList[iList]][sval];
					}
#					printf("%s\t%.2f\t%d\n", horMonList[iList], sval, monScoreA[horMonList[iList]][sval]) >> fbase "-score-stat.out";
					maxScore = (maxScore > sval) ? sval : maxScore;
					minScore = (minScore > sval) ? minScore : sval;
				}
				printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
					monScoreDistrA[horMonList[iList]]["leftCount"], \
					monScoreDistrA[horMonList[iList]]["middleCount"], \
					monScoreDistrA[horMonList[iList]]["rightCount"]) >> fbase "-score-stat.out";
			}
# Length distribution stats
			if (lengthstat != 0) {
				for (lval in monLengthA[horMonList[iList]]) {
					if (lval <= 50) {
						monLengthDistrA[horMonList[iList]]["leftCount"] += monLengthA[horMonList[iList]][lval];
					}
					if ((lval > 50)&&(lval < 76)) {
						monLengthDistrA[horMonList[iList]]["middleCount"] += monLengthA[horMonList[iList]][lval];
					}
					if (lval >= 76) {
						monLengthDistrA[horMonList[iList]]["rightCount"] += monLengthA[horMonList[iList]][lval];
					}
#					printf("%s\t%.2f\t%d\n", horMonList[iList], lval, monLengthA[horMonList[iList]][lval]) >> fbase "-length-stat.out";
					maxLength = (maxLength > lval) ? maxLength : lval;
					minLength = (minLength > lval) ? lval : minLength;
				}
				printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
					monLengthDistrA[horMonList[iList]]["leftCount"], \
					monLengthDistrA[horMonList[iList]]["middleCount"], \
					monLengthDistrA[horMonList[iList]]["rightCount"]) >> fbase "-length-stat.out";
			}
		} else {
			printf("%s\t%d\n", horMonList[iList], monA[horMonList[iList]]+0);
			if (scorestat != 0) {
				printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
					monScoreDistrA[horMonList[iList]]["leftCount"], \
					monScoreDistrA[horMonList[iList]]["middleCount"], \
					monScoreDistrA[horMonList[iList]]["rightCount"]) >> fbase "-score-stat.out";
			}
			if (lengthstat != 0) {
				printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
					monLengthDistrA[horMonList[iList]]["leftCount"], \
					monLengthDistrA[horMonList[iList]]["middleCount"], \
					monLengthDistrA[horMonList[iList]]["rightCount"]) >> fbase "-length-stat.out";
			}
		}
	}

# Print HOR groups
    print "\n#List_of_HORs_by_group";
    gNumsL = length(gNums);
    for (gNum=1; gNum<=gNumsL; gNum++) {
        if (gNum == 1) printf("%s\n","Group_1: Only \"live\" HORs");
        if (gNum == 2) printf("%s\n","Group_2: Only \"dead\" high copy HORs");
        if (gNum == 3) printf("%s\n","Group_3: Only \"dead\" mid copy HORs");
        if (gNum == 4) printf("%s\n","Group_4: Only \"dead\" low copy HORs");
        for (igh=1; igh<=gHorNum; igh++) {
            horComment = horComments[gNum][horGroups[gNum][igh]];
            if (horGroups[gNum][igh] in horList) {
                printf("%s%d_%s", "g", gNum, horGroups[gNum][igh]);
                if (horComment != "") printf("_%s", horComment);
                printf("\t%d\n", horA[horGroups[gNum][igh]]+0);
            }
        }
    }

# Print elements absent in HOR mon list
	for (iInput in monA) {
		if (!(iInput in horMonList)) {
			foundAbsentMons++;
			if (foundAbsentMons == 1) {
				print "\n#Elements_absent_in_HOR_monomers_list:";
			}
            absentMonA[iInput] = monA[iInput]+0;
#			printf("%s\t%d\n", iInput, monA[iInput]+0);
		}
	}
    absentN = asorti(absentMonA, absentMonAs);
    for (i=1; i<=absentN; i++) {
        printf("%s\t%d\n", absentMonAs[i], monA[absentMonAs[i]]+0);
    }

#	n = asorti(monA, monAs);
#	for (i=1; i<=n; i++) {
#		printf("%s\t%d\n", monAs[i], monA[monAs[i]]+0);
#		m += monA[monAs[i]];
#	}
	printf("%s\t%d\n", "Monomer_counter:", mCount+0);

	printf("%s\t%.2f\n", "Max_score:", maxScore);
	printf("%s\t%.2f\n", "Min_score:", minScore);
	printf("%s\t%d\n", "Unique_score_values:", length(scA));
	printf("%s\t%.2f\n", "Max_length:", maxLength);
	printf("%s\t%.2f\n", "Min_length:", minLength);
}
