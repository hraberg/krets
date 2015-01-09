# krets

Electronic circuit analysis in Clojure using Modified Nodal Analysis. Extremely rudimentary.

During 1996-97 I was writing a [Nord Lead clone in C++](http://www.student.nada.kth.se/~raberg/vl.html),  which kind of worked but was never properly finished. So almost 10 years later, I've decided to start working on another software synthesizer using a different approach. I'm learning a lot - I haven't looked at electronics since high school - and I'm mainly having fun exploring the subject.

The goal here is real-time circuit modelling of an analog synthesizer. A mildly realistic target is the [Music From Outer Space Noise Toaster](http://www.musicfromouterspace.com/analogsynth_new/NOISETOASTER/NOISETOASTER.php) which is the subject of Ray Wilson's book [Make: Analog Synthesizers](http://www.makershed.com/products/make-analog-synthesizers) A first step towards that is the [Alien Screamer](http://www.musicfromouterspace.com/analogsynth_new/ALIENSCREAMER/ALIENSCREAMER.php). This repository is very far from that goal though! I expect to ditch Clojure for some other language - Rust, Nim or C - once (and if) I figure out how it's supposed to work.

To achieve this we need to be able to simulate `V` voltage sources , `I` current sources, `E` voltage controlled voltage sources, `R` resistors, `C` capacitors, `U` op amps, `D` diodes, `Q` bipolar junction transistors and `J` junction gate field-effect transistors. Preferably in real-time. There's no intention (or knowledge) to build a full-fledged general purpose circuit simulator.

The industry standard circuit simulator is called Spice. 'Krets' is Swedish for circuit.

Tests and implementation guidance from http://www.ecircuitcenter.com/SPICEtopics.htm

## References

### Synthesizers

#### DIY Analogs

* http://musicfromouterspace.com
* http://rubidium.dyndns.org/~magnus/synths/friends/stopp/

#### Software Circuit Simulation

* http://www.u-he.com/cms/diva
* http://www.rolandus.com/blog/2014/02/14/analog-circuit-behavior-acb/
* http://www.arturia.com/cs-80v/tae%C2%AE
* http://www.cytomic.com/drop
* http://www.livespice.org/

Industry conference for plugin developers.

* http://www.dafx.de/

### Modified Node Analysis / Electronics

This is the technique used by Spice and most other simulators.

* http://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA_All.html
* http://www.ecircuitcenter.com/SPICEtopics.htm
* http://qucs.sourceforge.net/docs/technical/technical.pdf
* http://www.allaboutcircuits.com/
* https://ccrma.stanford.edu/~dtyeh/papers/DavidYehThesissinglesided.pdf

#### Spice

* http://bwrcs.eecs.berkeley.edu/Classes/IcBook/SPICE/
* http://www.linear.com/designtools/software/#LTspice
* http://newton.ex.ac.uk/teaching/CDHW/Electronics2/userguide/

#### Other MNA simulators

* http://www.falstad.com/circuit/ (Java)
* https://github.com/zupolgec/circuit-simulator (JavaScript)
* https://github.com/Qucs/qucs/ (C++)
* https://github.com/ahkab/ahkab (Python)

### Wave Digital Filters

Another interesting approach, needs more hands-on knowledge to manage circuits with multiple non-linear elements.

* https://ccrma.stanford.edu/~dtyeh/papers/wdftutorial.pdf
* http://www.nireaktor.com/reaktor-tutorials/wave-digital-filters-in-reaktor/
* http://home.deib.polimi.it/sarti/CV_and_publications_files/2010_IEEE_TrASLP_virtual_analog_modeling_WD.pdf
* http://legacy.spa.aalto.fi/software/BlockCompiler/

Yet another approach, based on random walks:

* http://www.ece.umn.edu/~sachin/conf/dac03.pdf

### VST Development

* http://www.juce.com/
* https://github.com/dwu/ClojureVST

### Dependencies

* https://github.com/mikera/core.matrix
* http://www.jfree.org/

## License

Copyright © 2015 Håkan Råberg

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.