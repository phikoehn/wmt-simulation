#!/usr/bin/perl -w

use strict;
use Math::Random;
use Getopt::Long "GetOptions";

my $STANDARD_DEVIATON = 10;
my $SYSTEM_COUNT = 15;

my @LOG_FACTORIAL = (0,0);
my @FACTORIAL = (1,1);
my $DEBUG = 0;
my $NON_RANDOM = 0;
my $COMPARE_RANKING_METHODS = 0;

my $MAX_CYCLE = 5;
my $USE_TARGETED_SAMPLING = 0; # 500
my $NUMBER_OF_EXPERIMENTS = 100;
my $SAMPLES = 5000;
my $STEP = 500;

my $REPORT_P_VALUE = 0;
my $REPORT_JUDGEMENT_COUNT = 0;
my $REPORT_CLUSTERS = 0;
my $REPORT_CLUSTERS_PAIRWISE = 0;
my $REPORT_RANKING = 0;

die unless &GetOptions('targeted-sampling=i' => \$USE_TARGETED_SAMPLING,
                       'runs=i' => \$NUMBER_OF_EXPERIMENTS,
		       'systems=i' => \$SYSTEM_COUNT,
                       'max-cycle=i' => \$MAX_CYCLE,
		       'non-random' => \$NON_RANDOM,
		       'samples=i' => \$SAMPLES,
		       'sd=f' => \$STANDARD_DEVIATON,
		       'compare-ranking-methods' => \$COMPARE_RANKING_METHODS,
                       'debug' => \$DEBUG,
		       'step=i' => \$STEP, 
                       'report-clusters' => \$REPORT_CLUSTERS,
                       'report-clusters-pairwise' => \$REPORT_CLUSTERS_PAIRWISE,
		       'report-ranking' => \$REPORT_RANKING,
                       'report-p-level' => \$REPORT_P_VALUE,
                       'report-judgement-count' => \$REPORT_JUDGEMENT_COUNT);

print "$SYSTEM_COUNT systems, average over $NUMBER_OF_EXPERIMENTS experiments with up to $SAMPLES samples, sd=$STANDARD_DEVIATON\n";

my @SYSTEM_MEDIAN;
my (%RANKING_QUALITY,%DISTINCTION_REPORT);
my (%CLUSTER,%VIOLATION,%RANGE_SIZE,%RANGE_VIOLATION);
my (%WIN,%TOTAL,@SAMPLE);
for(my $exp=0;$exp<$NUMBER_OF_EXPERIMENTS;$exp++) {
    (%WIN,%TOTAL,@SAMPLE,@SYSTEM_MEDIAN) = ((),(),(),());
    for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
	if ($NON_RANDOM) {
	    push @SYSTEM_MEDIAN, 10/$SYSTEM_COUNT*$i;
	}
	else {
	    push @SYSTEM_MEDIAN, rand(10);
	}
    }
    @SYSTEM_MEDIAN = reverse sort { $a <=> $b } @SYSTEM_MEDIAN;
    
    for(my $i=1;$i<=$SAMPLES;$i++) {
	if ($USE_TARGETED_SAMPLING>0 && $i > $USE_TARGETED_SAMPLING) {
	    &targeted_sample_judgment();
	}
	else {
	    &sample_judgment();
	}
	if ($i > 0 && $i % $STEP == 0) {
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
	    if ($COMPARE_RANKING_METHODS) {
		@RANKING = &rank_lopez_approximation();
		$RANKING_QUALITY{$i}{"lopez"} += &evaluate_ranking(@RANKING);
		#print "lopez: ".join(" ",@RANKING).": ".&evaluate_ranking(@RANKING)."\n" if $i == $SAMPLES;
	    }
	    
	    @RANKING = &rank_by_average_number_of_wins($i==$SAMPLES);
	    $RANKING_QUALITY{$i}{"exp_w"} += &evaluate_ranking(@RANKING);
	    #print "exp_w: ".join(" ",@RANKING).": ".&evaluate_ranking(@RANKING)."\n" if $i == $SAMPLES;

	    if ($COMPARE_RANKING_METHODS) {
		@RANKING = &rank_bojar();
		$RANKING_QUALITY{$i}{"bojar"} += &evaluate_ranking(@RANKING);
		#print "bojar: ".join(" ",@RANKING).": ".&evaluate_ranking(@RANKING)."\n" if $i == $SAMPLES;
	    }
		
	    if ($REPORT_CLUSTERS || $REPORT_CLUSTERS_PAIRWISE) {
		my ($range_size,$range_violation_count,$clusters,$violation_count) = $REPORT_CLUSTERS_PAIRWISE ? &cluster_based_on_significant_pairwise_distinctions() : &bootstrap_clusters();
		print "$i: $clusters clusters, $violation_count violation_count\n";
		$RANGE_SIZE{$i} += $range_size;
		$RANGE_VIOLATION{$i} += $range_violation_count;
		$CLUSTER{$i} += $clusters;
		$VIOLATION{$i} += $violation_count;
	    }
	}
    }
}

print "DISTINCTION_REPORT\n" if scalar keys %DISTINCTION_REPORT;
foreach my $i (sort { $a <=> $b } keys %DISTINCTION_REPORT) {
    printf "%5d ",$i;
    foreach my $p_level (sort { $a <=> $b } keys %{$DISTINCTION_REPORT{$i}}) {
	printf "%.02f: %4.1f\t",$p_level/100,$DISTINCTION_REPORT{$i}{$p_level}*100/$NUMBER_OF_EXPERIMENTS;
    }
    print "\n";
}

print "RANKING_QUALITY (number of pairwise flips)\n" if scalar keys %RANKING_QUALITY;
foreach my $i (sort { $a <=> $b } keys %RANKING_QUALITY) { 
	printf "%5d",$i;
	foreach my $method (keys %{$RANKING_QUALITY{$i}}) {
	    printf(" %s: %4.1f",$method,$RANKING_QUALITY{$i}{$method}/$NUMBER_OF_EXPERIMENTS);
	}
	print "\n";
}

print "CLUSTER".($REPORT_CLUSTERS_PAIRWISE?" (pairwise method)":"")."\n" if scalar keys %CLUSTER;
foreach my $i (sort { $a <=> $b } keys %CLUSTER) {
    printf "%5d: %.2f average range size, %.3f range violations, %.2f clusters, %.3f violations\n",$i,$RANGE_SIZE{$i}/$NUMBER_OF_EXPERIMENTS,$RANGE_VIOLATION{$i}/$NUMBER_OF_EXPERIMENTS/$SYSTEM_COUNT,$CLUSTER{$i}/$NUMBER_OF_EXPERIMENTS,$VIOLATION{$i}/$NUMBER_OF_EXPERIMENTS/$SYSTEM_COUNT;
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
    my ($final) = @_;
    my @SUM_WINS;
    for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
	for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
	    next if $i == $j;
	    next unless defined($WIN{$i}{$j});
	    $SUM_WINS[$i] += $WIN{$i}{$j} / $TOTAL{$i}{$j};
	}
    }
    my @RANKING = reverse sort { $SUM_WINS[$a] <=> $SUM_WINS[$b] } keys %TOTAL;
    if ($final && $REPORT_RANKING) {
	print "FINAL RANKING\n";
	foreach my $i (@RANKING) {
	    printf "%2d: %4.1f wins\n",$i,$SUM_WINS[$i]*100/$SYSTEM_COUNT;
	}
    }
    return @RANKING;
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
    &count_sample(\%SCORE,\%WIN,\%TOTAL);
    push @SAMPLE, \%SCORE;
}

sub count_sample {
    my ($SCORE,$WIN,$TOTAL) = @_;
    foreach my $system1 (keys %$SCORE) {
	foreach my $system2 (keys %$SCORE) {
	    next if $system1 eq $system2;
	    if ($$SCORE{$system1} > $$SCORE{$system2}) {
		$$WIN{$system1}{$system2}++;
	    }
	    $$TOTAL{$system1}{$system2}++;
	}
    }
}

sub resample {
    my ($WIN,$TOTAL) = @_;
    for(my $i=0;$i<scalar(@SAMPLE);$i++) {
	&count_sample($SAMPLE[rand(scalar(@SAMPLE))],$WIN,$TOTAL);
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
    &count_sample(\%SCORE,\%WIN,\%TOTAL);
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

sub cluster_based_on_significant_pairwise_distinctions {
    my (%RANGE);
    foreach my $system1 (sort keys %TOTAL) {
	my ($start,$end) = (0,$SYSTEM_COUNT-1);
	foreach my $system2 (sort keys %{$TOTAL{$system1}}) {
	    next if $system1 eq $system2;
	    my $win = 0;
	    $win = $WIN{$system1}{$system2} if defined($WIN{$system1}{$system2});
	    my $loss = $TOTAL{$system1}{$system2}-$win;
	    my $p_value = &compute_p_value($win,$loss);
	    if ($p_value<=0.05) {
		$start++ if $loss>$win;
		$end-- if $win>$loss;
	    }
	}
	$RANGE{$system1} = "$start $end";
    }
    return &common_cluster(\%RANGE);
}

sub bootstrap_clusters {
    # step 1: get rankings from resamples, count rank assignments per system
    my %RANK;
    my $num_resample = 1000;
    for(my $r=0;$r<$num_resample;$r++) {
	my (%WIN,%TOTAL);
	&resample(\%WIN,\%TOTAL);
	
	my @SUM_WINS;
	for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
	    for(my $j=0;$j<$SYSTEM_COUNT;$j++) {
		next if $i == $j;
		next unless defined($WIN{$i}{$j});
		$SUM_WINS[$i] += $WIN{$i}{$j} / $TOTAL{$i}{$j};
	    }
	}
	
	my @RANKING = reverse sort { $SUM_WINS[$a] <=> $SUM_WINS[$b] } keys %TOTAL;
	#print join(" ",@RANKING)."\n";
	for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
	    $RANK{$RANKING[$i]}{$i}++;
	}
    }
    
    # step 2: get rank ranges per system
    my %RANGE;
    for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
	my ($max_rank,$max_prob) = (-1,0);
	foreach my $rank (keys %{$RANK{$i}}) {
	    my $prob = $RANK{$i}{$rank};
	    if ($prob > $max_prob) {
		$max_prob = $prob;
		$max_rank = $rank;
	    }
	}
	my $total_prob = $max_prob;
	my ($start,$end) = ($max_rank,$max_rank);
	while($total_prob < $num_resample*.95) {
	    if (! defined($RANK{$i}{$start-1})) {
		$end++; 
		$total_prob += $RANK{$i}{$end};
	    }
	    elsif(! defined($RANK{$i}{$end+1})) {
		$start--;
		$total_prob += $RANK{$i}{$start};
	    }
	    elsif($RANK{$i}{$start-1} > $RANK{$i}{$end+1}) {
		$start--;
		$total_prob += $RANK{$i}{$start};          
	    }
	    else {
		$end++; 
		$total_prob += $RANK{$i}{$end};         
	    }
	}
	$RANGE{$i} = "$start $end";
    }
    return &common_cluster(\%RANGE);
}

sub common_cluster {
    my ($RANGE) = @_;

    my ($range_size,$range_violation_count,%INSEPERABLE_AFTER_RANK) = (0,0);
    for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
	my ($start,$end) = split(/ /,$$RANGE{$i});
	#print "$i: $start-$end\n";
	for(my $s=$start;$s<$end;$s++) {
	    $INSEPERABLE_AFTER_RANK{$s} = 1;
	}
	$range_size += ($end-$start)/$SYSTEM_COUNT;
	$range_violation_count++ if ($i < $start || $i > $end);
    }
    
    # step 3: assign to clusters
    my ($violation_count,%CLUSTER) = (0,());
    for(my $i=0;$i<$SYSTEM_COUNT;$i++) {
	my ($start,$end) = split(/ /,$$RANGE{$i});  
	while($start>0 && $INSEPERABLE_AFTER_RANK{$start-1}) {
	    $start--;
	}
	while($end<$SYSTEM_COUNT-1 && $INSEPERABLE_AFTER_RANK{$end}) {
	    $end++;
	}
	push @{$CLUSTER{$start}}, $i;
	$violation_count++ if $i < $start || $i > $end;
    }
    return ($range_size,$range_violation_count,scalar keys %CLUSTER,$violation_count);
    
    # step 4: output clusters
    foreach my $start (sort { $a <=> $b } keys %CLUSTER) {
	my $size = scalar @{$CLUSTER{$start}};
	printf ("%2d-%2d:",$start,$start + $size - 1);
	foreach my $i (@{$CLUSTER{$start}}) {
	    print " ".$i;
	    print "*" if ($i < $start || $i >= $start + $size);
	}
	print "\n";
    }     
}
