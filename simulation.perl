#!/usr/bin/perl -w

use strict;
use Math::Random;
use Getopt::Long "GetOptions";

my $STANDARD_DEVIATON = 10;
my $SYSTEM_COUNT = 15;

my @LOG_FACTORIAL = (0,0);
my @FACTORIAL = (1,1);
my $DEBUG = 0;

my $MAX_CYCLE = 5;
my $USE_TARGETED_SAMPLING = 0; # 500
my $NUMBER_OF_EXPERIMENTS = 1000;

my $REPORT_P_VALUE = 0;
my $REPORT_JUDGEMENT_COUNT = 0;

die unless &GetOptions('targeted-sampling=i' => \$USE_TARGETED_SAMPLING,
                       'runs=i' => \$NUMBER_OF_EXPERIMENTS,
                       'max-cycle=i' => \$MAX_CYCLE,
                       'debug' => \$DEBUG,
                       'report-p-level' => \$REPORT_P_VALUE,
                       'report-judgement-count' => \$REPORT_JUDGEMENT_COUNT);

my @SYSTEM_MEDIAN;
my (%RANKING_QUALITY,%DISTINCTION_REPORT);
my (%WIN,%TOTAL);
for(my $exp=0;$exp<$NUMBER_OF_EXPERIMENTS;$exp++) {
 (%WIN,%TOTAL,@SYSTEM_MEDIAN) = ((),(),());
 for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
  push @SYSTEM_MEDIAN, rand(10);
 }
 @SYSTEM_MEDIAN = reverse sort { $a <=> $b } @SYSTEM_MEDIAN;

 for(my $i=1;$i<=5000;$i++) {
  if ($USE_TARGETED_SAMPLING>0 && $i > $USE_TARGETED_SAMPLING) {
    &targeted_sample_judgment();
  }
  else {
    &sample_judgment();
  }
  if ($i > 0 && $i % 500 == 0) {
    my %DISTINCTION = &count_distinctions();
    foreach my $p_level (sort { $a <=> $b } keys %DISTINCTION) {
      $DISTINCTION_REPORT{$i}{$p_level} += $DISTINCTION{$p_level}/($SYSTEM_COUNT*($SYSTEM_COUNT-1));
    }

    if ($REPORT_P_VALUE || $REPORT_JUDGEMENT_COUNT) {
     print "   ";
     for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
      printf "%5d",$i;
     }
     print " ($i)\n";
     print "   ";
     for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
       printf " %4.2f", $SYSTEM_MEDIAN[$i];
     }
     print "\n";
     for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
      if ($REPORT_P_VALUE) {
        printf "%2d:",$i;
        for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
          if ($i == $j) { print "    -"; }
          else { printf " %4.2f",&compute_p_value($WIN{$i}{$j},$WIN{$j}{$i}); }
        }
        print "\n";
      }
      if ($REPORT_JUDGEMENT_COUNT) {
        printf "%2d:",$i;
        for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
          if ($i == $j) { print "    -"; }
          else { printf "%5d",$TOTAL{$i}{$j}; }
        }
        print "\n";
      }
     }
    }
    my @RANKING;
    @RANKING = &rank_lopez_approximation();
    $RANKING_QUALITY{$i}{"lopez"} += &evaluate_ranking(@RANKING);
    #print "lopez: ".join(" ",@RANKING).": ".&evaluate_ranking(@RANKING)."\n" if $i == 5000;

    @RANKING = &rank_by_average_number_of_wins();
    $RANKING_QUALITY{$i}{"exp_w"} += &evaluate_ranking(@RANKING);
    #print "exp_w: ".join(" ",@RANKING).": ".&evaluate_ranking(@RANKING)."\n" if $i == 5000;

    @RANKING = &rank_bojar();
    $RANKING_QUALITY{$i}{"bojar"} += &evaluate_ranking(@RANKING);
    #print "bojar: ".join(" ",@RANKING).": ".&evaluate_ranking(@RANKING)."\n" if $i == 5000;
  }
 }
}

foreach my $i (sort { $a <=> $b } keys %DISTINCTION_REPORT) {
  printf "%5d ",$i;
  foreach my $p_level (sort { $a <=> $b } keys %{$DISTINCTION_REPORT{$i}}) {
    printf "%.02f: %4.1f\t",$p_level/100,$DISTINCTION_REPORT{$i}{$p_level}*100/$NUMBER_OF_EXPERIMENTS;
  }
  print "\n";
}

foreach my $i (sort { $a <=> $b } keys %RANKING_QUALITY) { 
  printf "%5d",$i;
  foreach my $method (keys %{$RANKING_QUALITY{$i}}) {
    printf(" %s: %4.1f",$method,$RANKING_QUALITY{$i}{$method}/$NUMBER_OF_EXPERIMENTS);
  }
  print "\n";
}

sub rank_lopez_approximation {
  my %WEIGHT;
  for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
    for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
      $WIN{$j}{$i} = 0 unless defined($WIN{$j}{$i});
      $WIN{$i}{$j} = 0 unless defined($WIN{$i}{$j});
      if ($WIN{$j}{$i} > $WIN{$i}{$j}) {
        $WEIGHT{$i}{$j} = $WIN{$j}{$i} - $WIN{$i}{$j};
      }
    }
  }
  my @RANKING = &rank_bojar();
  print "bojar:        ".join(" ",@RANKING).", weight: ".&cost_ranking(\%WEIGHT,\@RANKING)."\n" if $DEBUG;
  &pairwise_hillclimbing(\%WEIGHT,\@RANKING);
  print "hillclimbing: ".join(" ",@RANKING).", weight: ".&cost_ranking(\%WEIGHT,\@RANKING)."\n" if $DEBUG;
  my @CYCLE;
  &detect_cycle(\%WEIGHT,\@RANKING,\@CYCLE);
  foreach (@CYCLE) {
    my ($min,$max) = split;
    if ($max-$min > $MAX_CYCLE) {
      print "CYCLE TOO BIG - APPROXIMATION\n" if $DEBUG;
      my $change = 1;
      while($change) {
        $change = 0;
        for(my $offset=0;$offset<=$max-$min-$MAX_CYCLE;$offset++) {
          print "(sub cycle ".($min+$offset)." ".($min+$offset+$MAX_CYCLE).")\n" if $DEBUG;
          $change++ if &resolve_cycle(\%WEIGHT,\@RANKING,$min+$offset,$min+$offset+$MAX_CYCLE);
        }
      }
    }
    else {
      print "(cycle $min $max)\n" if $DEBUG;
      &resolve_cycle(\%WEIGHT,\@RANKING,$min,$max);
    }
  }
  return @RANKING;
}

sub cost_ranking {
  my ($WEIGHT,$RANKING) = @_;
  my $cost = 0;

  for(my $i=0;$i<scalar(@$RANKING)-1;$i++) {
    for(my $j=$i+1;$j<scalar(@$RANKING)-1;$j++) {
      if (defined($$WEIGHT{$$RANKING[$i]}{$$RANKING[$j]})) {
        $cost += $$WEIGHT{$$RANKING[$i]}{$$RANKING[$j]};
      }
    }
  }
  return $cost;
}

sub pairwise_hillclimbing {
  my ($WEIGHT,$RANKING) = @_;
  my $change = 1;
  while($change) {
    $change = 0;
    for(my $i=0;$i<scalar(@$RANKING)-1;$i++) {
      if (defined($$WEIGHT{$$RANKING[$i]}{$$RANKING[$i+1]})) {
        my $tmp = $$RANKING[$i];
        $$RANKING[$i] = $$RANKING[$i+1];
        $$RANKING[$i+1] = $tmp;
        $change = 1;
	    }
    }
  }
}

sub detect_cycle {
  my ($WEIGHT,$RANKING,$CYCLE) = @_;
  for(my $i=0;$i<scalar(@$RANKING)-1;$i++) {
    for(my $j=$i+1;$j<scalar(@$RANKING)-1;$j++) {
      if (defined($$WEIGHT{$$RANKING[$i]}{$$RANKING[$j]})) {
        print "$i($$RANKING[$i]) worse than $j($$RANKING[$j]): ".$$WEIGHT{$$RANKING[$i]}{$$RANKING[$j]}."\n" if $DEBUG;
        my $overlapping = 0;
        foreach (@$CYCLE) {
          my ($from,$to) = split;
		      if (($from <= $i && $i <= $to) && ($from <= $j && $j <= $to)) {
			      $overlapping = 1;
		      }
		      elsif (($from <= $i && $i <= $to) || ($from <= $j && $j <= $to)) {
			      $from = $i if $i < $from;
			      $to = $j if $j > $to;
			      $_ = "$from $to";
			      $overlapping = 1;
		      }
		    }
		    push @$CYCLE,"$i $j" unless $overlapping;
#        print join(" ; ",@$CYCLE)."\n";
	    }
	  }
  }
}

sub resolve_cycle {
  my ($WEIGHT,$RANKING,$min,$max) = @_;
  my @BEST;  
  my $size = $max-$min+1;
  my $n_rankings = &factorial($size);
  my $best_cost = 9e9;
  for(my $r=0;$r<$n_rankings;$r++) {
	  my $rr = $r;
	  my @ASSIGNED;
	  my $cost = 0;
	  for(my $i=0;$i<$size;$i++) {
	    my $p;
	    {
		    use integer;
		    $p = $rr / &factorial($size-$i-1);
	    }
	    my $position = $p;
	    foreach (sort { $a <=> $b } @ASSIGNED) {
		    $position++ if $_ <= $position;
	    }
	    for(my $j=0;$j<$i;$j++) {
		    $cost += $$WEIGHT{$$RANKING[$ASSIGNED[$j]+$min]}{$$RANKING[$position+$min]} if defined($$WEIGHT{$$RANKING[$ASSIGNED[$j]+$min]}{$$RANKING[$position+$min]});
	    }
	    push @ASSIGNED,$position;
	    $rr -= $p * &factorial($size-$i-1);
	  }
    
    if ($r == 0) {
      $best_cost = $cost;
    }
	  elsif ($cost < $best_cost) {
	    print "$r:" if $DEBUG;
	    @BEST = ();
	    for(my $i=0;$i<$size;$i++) {
		    print " $$RANKING[$ASSIGNED[$i]+$min]" if $DEBUG;
		    push @BEST, $$RANKING[$ASSIGNED[$i]+$min];
	    }
	    print " $cost $best_cost\n" if $DEBUG;
	    $best_cost = $cost;
	  }
  }
  return 0 unless scalar @BEST;
  for(my $i=0;$i<$size;$i++) {
	  $$RANKING[$i+$min] = $BEST[$i];
  }
  print "corrected:       ".join(" ",@$RANKING).", weight: ".&cost_ranking($WEIGHT,$RANKING)."\n" if $DEBUG;
  return 1;
}


sub rank_lopez {
  my %WEIGHT;
  for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
    for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
      $WIN{$j}{$i} = 0 unless defined($WIN{$j}{$i});
      $WIN{$i}{$j} = 0 unless defined($WIN{$i}{$j});
      if ($WIN{$j}{$i} > $WIN{$i}{$j}) {
        $WEIGHT{$i}{$j} = $WIN{$j}{$i} - $WIN{$i}{$j};
      }
    }
  }
  my (%COST,%AGENDA,%BEST_FOR_SET);
  $COST{""} = 0;
  $AGENDA{""}++;
  my @RANKING;
  while(1) {
    my ($min_cost,$best) = (9e9);
    foreach my $ranking (keys %AGENDA) {
      if ($COST{$ranking} < $min_cost) {
        $min_cost = $COST{$ranking};
        $best = $ranking;
      }
    }
    delete $AGENDA{$best};

    print "$min_cost: $best\n";
    @RANKING = split(/ /,$best);
    last if scalar(@RANKING) == $SYSTEM_COUNT;

    my %ALREADY;
    foreach(@RANKING) {
      $ALREADY{$_}++;
    }

    for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
      next if defined($ALREADY{$i});
      my $cost = $min_cost;
      foreach my $old (@RANKING) {
        $cost += $WEIGHT{$old}{$i} if defined($WEIGHT{$old}) && defined($WEIGHT{$old}{$i});
      }
      my $new = "$best $i";
      my $set = join(" ",split(/ /,$new));
      next if defined($BEST_FOR_SET{$set}) && $cost > $COST{$BEST_FOR_SET{$set}};
      if (defined($BEST_FOR_SET{$set})) {
        delete $AGENDA{$BEST_FOR_SET{$set}} if defined($AGENDA{$BEST_FOR_SET{$set}});
        delete $COST{$set} if defined $COST{$set};
      }
      $COST{$new} = $cost;
      $BEST_FOR_SET{$set} = $new;
      $AGENDA{$new}++;;
      #print "$new: $cost\n";
    }
  }
  return @RANKING;
}

sub rank_by_average_number_of_wins {
  my @SUM_WINS;
  for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
    for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
      next if $i == $j;
      next unless defined($WIN{$i}{$j});
      $SUM_WINS[$i] += $WIN{$i}{$j} / $TOTAL{$i}{$j};
    }
  }
  return reverse sort { $SUM_WINS[$a] <=> $SUM_WINS[$b] } keys %TOTAL;
}

sub rank_bojar {
  my @RATIO;
  for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
    my ($wins,$total) = (0,0);
    for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
      next if $i == $j;
      $total += $TOTAL{$i}{$j};
      next unless defined($WIN{$i}{$j});
      $wins += $WIN{$i}{$j};
    }
    $RATIO[$i] = $wins/$total;
  }
  return reverse sort { $RATIO[$a] <=> $RATIO[$b] } keys %TOTAL;
}

sub evaluate_ranking{
  my $mismatch = 0;
  for(my $i=0; $i<$SYSTEM_COUNT; $i++) {
    $mismatch += abs($i-$_[$i]);
  }
  return $mismatch;
}

sub sample_judgment {
  my %SCORE;
  for(my $i=0;$i<5;$i++) {
    my $system;
    while(!defined($system) || defined($SCORE{$system})) {
      $system = int(rand($SYSTEM_COUNT));
    }
    $SCORE{$system} = random_normal(1, $SYSTEM_MEDIAN[$system], $STANDARD_DEVIATON);
  }
  &count_sample(\%SCORE);
}

sub count_sample {
  my ($SCORE) = @_;
  foreach my $system1 (keys %$SCORE) {
    foreach my $system2 (keys %$SCORE) {
      next if $system1 eq $system2;
      if ($$SCORE{$system1} > $$SCORE{$system2}) {
        $WIN{$system1}{$system2}++;
      }
      $TOTAL{$system1}{$system2}++;
    }
  }
}

sub targeted_sample_judgment {
  my %SCORE;
  my $system = int(rand($SYSTEM_COUNT));
  $SCORE{$system} = random_normal(1, $SYSTEM_MEDIAN[$system], $STANDARD_DEVIATON);
  for(my $i=0;$i<4;$i++) {
    my %CANDIDATE;
    foreach my $old (keys %SCORE) {
      for(my $s=0;$s<$SYSTEM_COUNT;$s++) {
        if (!defined($SCORE{$s}) && &compute_p_value($WIN{$old}{$s},$WIN{$s}{$old}) > 0.01) {
          $CANDIDATE{$s}++;
        }
      }
    }
    if (scalar(keys %CANDIDATE) > 0) {
      my @C = keys %CANDIDATE;
      my $system = $C[rand(scalar(@C))];
      $SCORE{$system} = random_normal(1, $SYSTEM_MEDIAN[$system], $STANDARD_DEVIATON);
    }
    else {
      my $system;
      while(!defined($system) || defined($SCORE{$system})) {
        $system = int(rand($SYSTEM_COUNT));
      }
      $SCORE{$system} = random_normal(1, $SYSTEM_MEDIAN[$system], $STANDARD_DEVIATON);
    }
  }
  &count_sample(\%SCORE);
}

sub compute_p_value {
  my ($win,$loss) = @_;
  my $total = $win + $loss;
  my $min = $win < $loss ? $win : $loss;
  my $p = 0;
  for (my $i=0;$i<=$min;$i++) {
    $p += exp(&binomial($i,$total));
  }
  return $p;
}

sub binomial {
  my ($k,$n) = @_;
  my $log_prob = &log_factorial($n) - &log_factorial($n-$k) - &log_factorial($k);
  $log_prob += $n * log(0.5);
#  print "($k,$n) = ".exp($log_prob)."\n";
  return $log_prob;
}

sub log_factorial {
  my ($x) = @_;
  if ($x > $#LOG_FACTORIAL) {
    for(my $i=$#LOG_FACTORIAL;$i<=$x;$i++) {
      $LOG_FACTORIAL[$i] = $LOG_FACTORIAL[$i-1] + log($i);
    }
  }
  return $LOG_FACTORIAL[$x];
}

sub factorial {
  my ($x) = @_;
  if ($x > $#FACTORIAL) {
    for(my $i=$#FACTORIAL;$i<=$x;$i++) {
      $FACTORIAL[$i] = $FACTORIAL[$i-1] * $i;
    }
  }
  return $FACTORIAL[$x];
}

sub count_distinctions {
  my %DISTINCTION;
  foreach my $system1 (sort keys %TOTAL) {
    foreach my $system2 (sort keys %{$TOTAL{$system1}}) {
      next if $system1 eq $system2;
      my $win = 0;
      $win = $WIN{$system1}{$system2} if defined($WIN{$system1}{$system2});
      my $loss = $TOTAL{$system1}{$system2}-$win;
      my $p_value = &compute_p_value($win,$loss);
      #print "$task $system1 $system2 $p_value\n" if $ratio == 1;
      if($p_value <= 0.1) {
        $DISTINCTION{10}++;
        if($p_value <= 0.05) {
          $DISTINCTION{5}++;
          if ($p_value <= 0.01) {
            $DISTINCTION{1}++;
          }
        }
      }
    }
  }
  return %DISTINCTION;
}
