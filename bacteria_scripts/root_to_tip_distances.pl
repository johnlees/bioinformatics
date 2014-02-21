#!/usr/bin/perl -w

#****************************************************************************************#
#* root_to_tip_distances.pl                                                             *#
#* Gets the distance along branches when travelling between the root (tip of the rooted *#
#* node) and each of the tips                                                           *#
#* Root name is hardcoded for now (M66), so specific to Cholera                         *#
#****************************************************************************************#

use strict;

use Getopt::Long;
use Bio::TreeIO;

#* 
#* Main
#*

#* Get input tree, and read into Bio::TreeIO data structure
my $input_file;
GetOptions ("input|i=s"  => \$input_file,
		   ) or die($!);

my $input_tree = new Bio::TreeIO(-file   => $input_file,
                                 -format => "nexus");
my $tree = $input_tree->next_tree;

#* This reroots a tree, based on the node name
#*
# $tree->reroot($tree->find_node(-id => 'M66'));

#* This prints the tree as an SVG
#*
# my $out = new Bio::TreeIO(-file => '>mytree.svg',
#                           -format => 'svggraph');
# $out->write_tree($tree);

#* Get the name of the root node
my $root = $tree->get_root_node;
print ("root=" . $root->id . "\n");

my $M66 = $tree->find_node(-id => 'M66');

#* Cycle through nodes, only processing those which are tips (based on whether they are
#* named or not
for my $node ( $tree->get_nodes )
{
	if (defined($node->id))
	{
		#* Remove the quotation marks around the node IDs
		my $node_name;
		if ($node->id =~ /'(.+)'/)
		{
			$node_name = $1;
		}
		else
		{
			$node_name = $node->id;
		}
		
		#* Get the distance between the tip and M66
		my $distance = $tree->distance(-nodes => [$node, $M66]);
		print $node_name . " $distance\n";
	}
}

exit(0);
