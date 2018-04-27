# Collect stats from BED files after file(s) has been passed through hmmertblout2bed.awk
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
# Usage: awk -v l="AS-SF-HORs-list.txt" -v scorestat=1 -v lengthstat=1 -f bedhorstats.awk input.bed > output.stats

BEGIN {
# HOR list preparation
  if(l) {
    while((getline __input < l) > 0) {
      # add all names from the list
      if(__input !~ /^g|^#/) {
        fullName = __input;
        if(fullName ~ /\./) {
          dotpos        = index(fullName, ".");
          nameBeforeDot = substr(fullName, 1, dotpos-1)
        } else {
          nameBeforeDot = fullName
        }
        if(!(nameBeforeDot in horList)) {
          horList[++hNum]        = nameBeforeDot;
          horList[nameBeforeDot] = hNum;
        }
        horMonList[++mNum]   = fullName;
        horMonList[fullName] = mNum
      }
      # add groups from the list
      if(__input ~ /^g/) {
        fullName = __input;
        split(fullName, tmpA);
        gNum = substr(tmpA[1], 2);
        gNums[gNum]++;
        gHorName   = tmpA[2];
        horComment = tmpA[3];
        horGroups[gNum][++gHorNum]  = gHorName;
        horGroups[gNum][gHorName]   = gHorNum;
        horComments[gNum][gHorName] = horComment
      }
    }
    close(l)
  }
}

NR==1 {
  # get file path, name and extention
  fpath = fbase = fext = FILENAME;
  sub(/\/[^/]*$/,  "", fpath);
  sub(/.*[/]/,     "", fbase);
  sub(/[.][^.]+$/, "", fbase);
  sub(/.*[.]/,     "", fext );

  if(fpath == FILENAME) fpath = "."
}

{
  # get data from BED
  fullName = $4;
  if(fullName ~ /\./) {
    dotpos        = index(fullName, ".");
    nameBeforeDot = substr(fullName, 1, dotpos-1)
  } else {
    nameBeforeDot = fullName
  }
  start  = $2;
  finish = $3;
  score  = $5;
  monLength = finish - start;

  horA[nameBeforeDot]++;
  monA[fullName]++;
  horScoreA[nameBeforeDot][score]++;
  horLengthA[nameBeforeDot][monLength]++;
  monScoreA[fullName][score]++;
  monLengthA[fullName][monLength]++;

  scA[sprintf("%d", score)];

  if(!(once++)) {
    maxScore  = minScore  = score;
    maxLength = minLength = monLength;
  }
}

END {
  fsend = "-score-stat.out";
  flend = "-length-stat.out";

  # print elements in SF
  print "#" fbase;
  print "#Monomers_by_Suprachomosomal_Family";
  if(scorestat) {
    printf("#%s\t\t\n", fbase) > fpath "/" fbase fsend;
    printf("%s\n", "#Distribution_of_scores") >> fpath "/" fbase fsend;
    printf("%s\t%s\t%s\t%s\n", "HOR_list", "S<=61", "61<S<97", "S>=97") >> \
           fpath "/" fbase fsend;
  }
  if(lengthstat) {
    printf("#%s\t\t\n", fbase) > fpath "/" fbase flend;
    printf("%s\n", "#Distribution_of_lengths") >> fpath "/" fbase flend;
    printf("%s\t%s\t%s\t%s\n", "HOR_list", "L<=50", "50<L<76", "L>=76") >> \
           fpath "/" fbase flend;
  }
  for(iInput in horA) {
    if(iInput ~ /^S1/) sf1 += horA[iInput];
    if(iInput ~ /^S2/) sf2 += horA[iInput];
    if(iInput ~ /^S3/) sf3 += horA[iInput];
    if(iInput ~ /^S5/) sf5 += horA[iInput];
    if(iInput ~ /^S6/) sf6 += horA[iInput];
#   if(iInput ~ /^S4\/6/) sf4_6 += horA[iInput];
    if(iInput ~ /^S4/) sf4 += horA[iInput];
  }
  printf("%s\t%d\n", "SF1", sf1+0);
  printf("%s\t%d\n", "SF2", sf2+0);
  printf("%s\t%d\n", "SF3", sf3+0);
  printf("%s\t%d\n", "SF5", sf5+0);
  printf("%s\t%d\n", "SF4", sf4+0);
  printf("%s\t%d\n", "SF6", sf6+0);

  # print elements in HOR list
  print "#HOR_list";
  for(iList=1; iList<=hNum; iList++) {
    if(horList[iList] in horA) {
      if(horA[horList[iList]]+0 > 0) {
        hCount++;
      }
      printf("%s\t%d\n", horList[iList], horA[horList[iList]]+0);
      # print score distribution
      if(scorestat) {
        for(sval in horScoreA[horList[iList]]) {
          if(sval <= 61) {
            horScoreDistrA[horList[iList]]["leftCount"] += \
                          horScoreA[horList[iList]][sval];
            totalHorScoreA["leftCount"] += horScoreA[horList[iList]][sval];
          }
          if((sval > 61)&&(sval < 97)) {
            horScoreDistrA[horList[iList]]["middleCount"] += \
                          horScoreA[horList[iList]][sval];
            totalHorScoreA["middleCount"] += horScoreA[horList[iList]][sval];
          }
          if(sval >= 97) {
            horScoreDistrA[horList[iList]]["rightCount"] += \
                          horScoreA[horList[iList]][sval];
            totalHorScoreA["rightCount"] += horScoreA[horList[iList]][sval];
          }
          if(debug) {
            printf("%s\t%.2f\t%d\n", horList[iList], sval, \
            horScoreA[horList[iList]][sval]) >> fpath "/" fbase fsend;
          }
        }
        printf("%s\t%d\t%d\t%d\n", horList[iList], \
          horScoreDistrA[horList[iList]]["leftCount"], \
          horScoreDistrA[horList[iList]]["middleCount"], \
          horScoreDistrA[horList[iList]]["rightCount"]) >> \
          fpath "/" fbase fsend;
      }
      # print length distribution
      if(lengthstat) {
        for(lval in horLengthA[horList[iList]]) {
          if(lval <= 50) {
            horLengthDistrA[horList[iList]]["leftCount"] += \
                           horLengthA[horList[iList]][lval];
            totalHorLengthA["leftCount"] += horLengthA[horList[iList]][lval];
          }
          if((lval > 50)&&(lval < 76)) {
            horLengthDistrA[horList[iList]]["middleCount"] += \
                           horLengthA[horList[iList]][lval];
            totalHorLengthA["middleCount"] += horLengthA[horList[iList]][lval];
          }
          if(lval >= 76) {
            horLengthDistrA[horList[iList]]["rightCount"] += \
                           horLengthA[horList[iList]][lval];
            totalHorLengthA["rightCount"] += horLengthA[horList[iList]][lval];
          }
          if(debug) {
            printf("%s\t%.2f\t%d\n", horList[iList], lval, \
            horLengthA[horList[iList]][lval]) >> fpath "/" fbase flend;
          }
        }
        printf("%s\t%d\t%d\t%d\n", horList[iList], \
          horLengthDistrA[horList[iList]]["leftCount"], \
          horLengthDistrA[horList[iList]]["middleCount"], \
          horLengthDistrA[horList[iList]]["rightCount"]) >> \
          fpath "/" fbase flend;
      }
    } else {
      printf("%s\t%d\n", horList[iList], horA[horList[iList]]+0);
      if(scorestat) {
        printf("%s\t%d\t%d\t%d\n", horList[iList], \
          horScoreDistrA[horList[iList]]["leftCount"], \
          horScoreDistrA[horList[iList]]["middleCount"], \
          horScoreDistrA[horList[iList]]["rightCount"]) >> \
          fpath "/" fbase fsend;
      }
      if(lengthstat) {
        printf("%s\t%d\t%d\t%d\n", horList[iList], \
          horLengthDistrA[horList[iList]]["leftCount"], \
          horLengthDistrA[horList[iList]]["middleCount"], \
          horLengthDistrA[horList[iList]]["rightCount"]) >> \
          fpath "/" fbase flend;
      }
    }
  }

  # print total score values
  if(scorestat) {
    printf("%s\t%d\t%d\t%d\n", "Total_in_group:", \
      totalHorScoreA["leftCount"], \
      totalHorScoreA["middleCount"], \
      totalHorScoreA["rightCount"]) >> fpath "/" fbase fsend;
  }
  # print total length values
  if(lengthstat) {
    printf("%s\t%d\t%d\t%d\n", "Total_in_groups:", \
      totalHorLengthA["leftCount"], \
      totalHorLengthA["middleCount"], \
      totalHorLengthA["rightCount"]) >> fpath "/" fbase flend;
  }
  # print items absent in HOR list
  for(iInput in horA) {
    if(!(iInput in horList)) {
      foundAbsentHORs++;
      if(foundAbsentHORs == 1) {
        print "\n#Items_absent_in_HOR_list:";
      }
      printf("%s\t%d\n", iInput, horA[iInput]+0);
    }
  }

  printf("%s\t%d\n\n", "HORs_with_non-zero_amount_of_monomers:", hCount+0);

  print "#List_of_HOR_monomers";
  if(scorestat) {
    printf("\n%s\t%s\t%s\t%s\n", \
      "List_of_HOR_monomers", \
      "S<=61", "61<S<97", "S>=97") >> fpath "/" fbase fsend;
  }
  if(lengthstat) {
    printf("\n%s\t%s\t%s\t%s\n", \
      "List_of_HOR_monomers", \
      "L<=50", "50<L<76", "L>=76") >> fpath "/" fbase flend;
  }

  # print elements in HOR monomer list
  for(iList=1; iList<=mNum; iList++) {
    if(horMonList[iList] in monA) {
      mCount += monA[horMonList[iList]];
      printf("%s\t%d\n", horMonList[iList], monA[horMonList[iList]]+0);
      # score distribution stats
      if(scorestat) {
        for(sval in monScoreA[horMonList[iList]]) {
          if(sval <= 61) {
            monScoreDistrA[horMonList[iList]]["leftCount"] += \
                          monScoreA[horMonList[iList]][sval];
          }
          if((sval > 61)&&(sval < 97)) {
            monScoreDistrA[horMonList[iList]]["middleCount"] += \
                          monScoreA[horMonList[iList]][sval];
          }
          if(sval >= 97) {
            monScoreDistrA[horMonList[iList]]["rightCount"] += \
                          monScoreA[horMonList[iList]][sval];
          }
          if(debug) {
            printf("%s\t%.2f\t%d\n", horMonList[iList], sval, \
            monScoreA[horMonList[iList]][sval]) >> fpath "/" fbase fsend;
          }
          # get min/max score values
          maxScore = (sval+0 > maxScore+0) ? sval : maxScore;
          minScore = (sval+0 < minScore+0) ? sval : minScore;
        }
        printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
          monScoreDistrA[horMonList[iList]]["leftCount"], \
          monScoreDistrA[horMonList[iList]]["middleCount"], \
          monScoreDistrA[horMonList[iList]]["rightCount"]) >> \
          fpath "/" fbase fsend;
      }
      # length distribution stats
      if(lengthstat) {
        for(lval in monLengthA[horMonList[iList]]) {
          if(lval <= 50) {
            monLengthDistrA[horMonList[iList]]["leftCount"] += \
                           monLengthA[horMonList[iList]][lval];
          }
          if((lval > 50)&&(lval < 76)) {
            monLengthDistrA[horMonList[iList]]["middleCount"] += \
                           monLengthA[horMonList[iList]][lval];
          }
          if(lval >= 76) {
            monLengthDistrA[horMonList[iList]]["rightCount"] += \
                           monLengthA[horMonList[iList]][lval];
          }
          if(debug) {
            printf("%s\t%.2f\t%d\n", horMonList[iList], lval, \
              monLengthA[horMonList[iList]][lval]) >> fpath "/" fbase flend;
          }
          # get min/max length values
          maxLength = (lval+0 > maxLength+0) ? lval : maxLength;
          minLength = (lval+0 < minLength+0) ? lval : minLength;
        }
        printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
          monLengthDistrA[horMonList[iList]]["leftCount"], \
          monLengthDistrA[horMonList[iList]]["middleCount"], \
          monLengthDistrA[horMonList[iList]]["rightCount"]) >> \
          fpath "/" fbase flend;
      }
    } else {
      printf("%s\t%d\n", horMonList[iList], monA[horMonList[iList]]+0);
      if(scorestat) {
        printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
          monScoreDistrA[horMonList[iList]]["leftCount"], \
          monScoreDistrA[horMonList[iList]]["middleCount"], \
          monScoreDistrA[horMonList[iList]]["rightCount"]) >> \
          fpath "/" fbase fsend;
      }
      if(lengthstat) {
        printf("%s\t%d\t%d\t%d\n", horMonList[iList], \
          monLengthDistrA[horMonList[iList]]["leftCount"], \
          monLengthDistrA[horMonList[iList]]["middleCount"], \
          monLengthDistrA[horMonList[iList]]["rightCount"]) >> \
          fpath "/" fbase flend;
      }
    }
  }

  # print HOR groups
  print "\n#HOR_groups";
  gNumsL = length(gNums);
  for(gNum=1; gNum<=gNumsL; gNum++) {
    if(gNum == 1) printf("%s\n","Group#1: live HORs");
    if(gNum == 2) printf("%s\n","Group#2: high copy HORs");
    if(gNum == 3) printf("%s\n","Group#3: average copy HORs");
    if(gNum == 4) printf("%s\n","Group#4: low copy HORs");
    for(igh=1; igh<=gHorNum; igh++) {
      horComment = horComments[gNum][horGroups[gNum][igh]];
      if(horGroups[gNum][igh] in horList) {
        printf("%s%d_%s", "g", gNum, horGroups[gNum][igh]);
        if(horComment != "") printf("_%s", horComment);
          printf("\t%d\n", horA[horGroups[gNum][igh]]+0);
      }
    }
  }

  # print items absent in HOR mon list
  for(iInput in monA) {
    if(!(iInput in horMonList)) {
      foundAbsentMons++;
      if(foundAbsentMons == 1) {
        print "\n#Items_absent_in_HOR_monomers_list:";
      }
      absentMonA[iInput] = monA[iInput]+0;
      if(debug) {
        printf("%s\t%d\n", iInput, monA[iInput]+0);
      }
    }
  }

  absentN = asorti(absentMonA, absentMonAs);
  for(i=1; i<=absentN; i++) {
    printf("%s\t%d\n", absentMonAs[i], monA[absentMonAs[i]]+0);
  }

  printf("%s\t%d\n",   "Monomer_counter:",     mCount+0);
  printf("%s\t%.2f\n", "Min_score:",           minScore);
  printf("%s\t%.2f\n", "Max_score:",           maxScore);
  printf("%s\t%d\n",   "Unique_score_values:", length(scA));
  printf("%s\t%.2f\n", "Min_length:",          minLength);
  printf("%s\t%.2f\n", "Max_length:",          maxLength);
}
