#! perl -w

use strict;
use warnings;

my $funclist = "funclist.txt";
my $tsvfile = "funclist.tsv";
my $graphmlFile = "funclist.graphml";

my %colHash = (
    '00-generics'   => "#85660D", # goldenrod4
    '01-cloneMaker' => "#1CBE4F", # springreen3
    '02-prefit'     => "#90AD1C", # darkolivegreen3
    '03-psi'        => "#16FF32", # green
    algorithm     => "#325A9B", # dodgerblue4
    functions     => "#2ED9FF", # lightskyblue
    running       => "#822E1C", # tomato4
    simgen        => "#F6222E", # firebrick1
    assessing     => "#683B79", # mediumorchid4
    clinical      => "#C4451C", # tomato3
    cnv_survival  => "#FEAF16", # darkgoldenrod1
    figures       => "#1C8356", # seagreen
    );

my %fgHash = (
    '00-generics'  => "#FFFFFF",
    '01-cloneMaker'=> "#000000",
    '02-prefit'    => "#000000",
    '03-psi'       => "#000000",
    algorithm      => "#FFFFFF",
    functions      => "#000000",
    running        => "#FFFFFF",
    simgen         => "#000000",
    assessing      => "#FFFFFF",
    clinical       => "#FFFFFF",
    cnv_survival   => "#000000",
    figures        => "#FFFFFF",
    );

## GraphML outer structure
my $gmlHead = <<'ENDHEAD'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:java="http://www.yworks.com/xml/yfiles-common/1.0/java" xmlns:sys="http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0" xmlns:x="http://www.yworks.com/xml/yfiles-common/markup/2.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:y="http://www.yworks.com/xml/graphml" xmlns:yed="http://www.yworks.com/xml/yed/3" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd">
  <!--Created by yEd 3.15.0.2-->
  <key attr.name="Description" attr.type="string" for="graph" id="d0"/>
  <key for="port" id="d1" yfiles.type="portgraphics"/>
  <key for="port" id="d2" yfiles.type="portgeometry"/>
  <key for="port" id="d3" yfiles.type="portuserdata"/>
  <key attr.name="url" attr.type="string" for="node" id="d4"/>
  <key attr.name="description" attr.type="string" for="node" id="d5"/>
  <key for="node" id="d6" yfiles.type="nodegraphics"/>
  <key for="graphml" id="d7" yfiles.type="resources"/>
  <key attr.name="url" attr.type="string" for="edge" id="d8"/>
  <key attr.name="description" attr.type="string" for="edge" id="d9"/>
  <key for="edge" id="d10" yfiles.type="edgegraphics"/>
  <graph edgedefault="directed" id="G">
    <data key="d0"/>
ENDHEAD
    ;
my $gmlTail = <<'ENDTAIL'
  </graph>
  <data key="d7">
    <y:Resources/>
  </data>
</graphml>
ENDTAIL
    ;


## Node Structure
my $beginNode =<<'BNODE'
      <data key="d5"/>
      <data key="d6">
        <y:ShapeNode>
          <y:Geometry height="30.0" width="100.033203125" x="467.9833984375" y="91.0"/>
          <y:Fill color=
BNODE
    ;
my $midNode =<<'MNODE'
 transparent="false"/>
          <y:BorderStyle color="#000000" type="line" width="1.0"/>
          <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="18.701171875" modelName="custom" textColor=
MNODE
    ;
my $mid2node =<<'M2'
 visible="true" width="90.033203125" x="5.0" y="5.6494140625">
M2
    ;
my $endNode =<<'ENODE'
<y:LabelModel>
              <y:SmartNodeLabelModel distance="4.0"/>
            </y:LabelModel>
            <y:ModelParameter>
              <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
            </y:ModelParameter>
          </y:NodeLabel>
          <y:Shape type="roundrectangle"/>
        </y:ShapeNode>
      </data>
    </node>
ENODE
    ;

## Edge Structure
#    <edge id="e0" source="n112" target="n47">
my $makeEdge = <<'MEDGE'
      <data key="d9"/>
      <data key="d10">
        <y:PolyLineEdge>
          <y:Path sx="0.0" sy="0.0" tx="0.0" ty="0.0"/>
          <y:LineStyle color="#000000" type="line" width="1.0"/>
          <y:Arrows source="none" target="standard"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>
MEDGE
    ;


my %nodeHash = ();

open(GML, ">$graphmlFile") or die "Unable to create '$graphmlFile': $!\n";
print GML $gmlHead;

open(FUNC, "<$funclist") or die "Unable to open '$funclist': $!\n";
open(TSV, ">$tsvfile") or die "Unable to create '$tsvfile': $!\n";
my $nodeCount = 0;
while(my $func = <FUNC>) {
    chomp($func);
    next unless($func);    # skip empty lines
    next if($func =~ /^#/); # skip over comments
    if ($func =~ /^(.*?)\.R\:(.*?)<-/) { # filename.R:[optional white space]funcname[opt ws]<-
	my $filename = $1;
	my $funcname = $2;
	$funcname =~ s/\s//g;
	print TSV "$filename\t$funcname\t$nodeCount\n";
	my $nodeNumber = "    <node id=\"n$nodeCount\">\n";
	$nodeHash{$funcname} = "n$nodeCount";
	++$nodeCount;
	my $nodeColor = $colHash{$filename} || "#999999";
	my $txtColor = $fgHash{$filename} || "#FF0000";
	print GML $nodeNumber, $beginNode, "\"$nodeColor\"", 
                  $midNode, "\"$txtColor\"", $mid2node,
	          "$filename: $funcname", $endNode;
    } else {
	print STDERR "Bad line: $func\n";
    }
}
close(TSV);
close(FUNC);

if (0) {
    print STDERR "Node Names:\n";
    foreach my $key (sort keys %nodeHash) {
	print STDERR "'$key' => $nodeHash{$key}\n";
    }
}

#    <edge id="e0" source="n112" target="n47">
my $edgeset = "edgeSet.txt";
open(EDGE, "<$edgeset") or die "Unable to open '$edgeset': $!\n";
my $edgeCount = 0;
while(my $edge = <EDGE>) {
    chomp $edge;
    my ($begin, $end) = split /\|/, $edge;
    my $caller = $nodeHash{$begin};
#    print "Trying to build edge from '$caller'\n";
    die "Bad begin funcname: '$begin'\n" unless($caller);
    my $callee = $nodeHash{$end} or die "Bad end funcname: '$end'\n";
    my $line = "   <edge id =\"e$edgeCount\" source=\"$caller\" target=\"$callee\">\n";
    ++$edgeCount;
    print GML $line, $makeEdge;
}
close(EDGE);

print GML $gmlTail;
close(GML);

exit;
__END__
