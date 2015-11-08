# Bionode-seq
A module for DNA, RNA and protein sequences manipulation.

## Install and load
This is a very simple module and doesn't have a command line interface yet as it's mostly used in the browser.
You can load and install it like this:

```bash
npm install bionode-seq
```

```javascript
var seq = require('bionode-seq')
```

### Usage

### Check sequence type

Takes a sequence string and checks if it's DNA, RNA or protein (returns 'dna', 'rna', 'protein' or undefined). Other optional arguments include threshold, length and index (see below).
```javascript
 seq.checkType("ATGACCCTGAGAAGAGCACCG");
 //=> "dna"
 seq.checkType("AUGACCCUGAAGGUGAAUGAA");
 //=> "rna"
 seq.checkType("MAYKSGKRPTFFEVFKAHCSDS");
 //=> "protein"
 seq.checkType("1234567891234567ATGACC");
 //=> undefined
```
 By default, the method has a 90% threshold, however, this can be altered as required.

     seq.checkType("1234567891234567ATGACC", 0.8);
     => undefined
     seq.checkType("--------MAYKSGKRPTFFEV", 0.7);
     => "protein"

 The length value specifies the length of the sequence to be analyse (default 10000). If your sequence is extremely long, you may want to analyse a shorter sub-section to reduce the computational burden.

     seq.checkType('A Very Long Sequence', 0.9, 1000);
     => Type based on the first 1000 characters

 The index value specifies the point on the sequence from which the sequence is to be analysed. Perhaps you know that there are lot of gaps at the start of the sequence.

     seq.checkType("--------MAYKSGKRPTFFEV", 0.9, 10000, 8);
     => "protein"

 Takes a sequence type argument and returns a function to complement bases.
 ### Reverse sequence
 Takes sequence string and returns the reverse sequence.

     seq.reverse("ATGACCCTGAAGGTGAA");
     => "AAGTGGAAGTCCCAGTA"
 ### (Reverse) complement sequence
 Takes a sequence string and optional boolean for reverse, and returns its complement.

     seq.complement("ATGACCCTGAAGGTGAA");
     => "TACTGGGACTTCCACTT"
     seq.complement("ATGACCCTGAAGGTGAA", true);
     => "TTCACCTTCAGGGTCAT"
     Alias
     seq.reverseComplement("ATGACCCTGAAGGTGAA");
     => "TTCACCTTCAGGGTCAT"
 Takes a sequence string and returns the reverse complement (syntax sugar).
 ### Transcribe base
 Takes a base character and returns the transcript base.

     seq.getTranscribedBase("A");
     => "U"
     seq.getTranscribedBase("T");
     => "A"
     seq.getTranscribedBase("t");
     => "a"
     seq.getTranscribedBase("C");
     => "G"
 ### Get codon amino acid
 Takes an RNA codon and returns the translated amino acid.

     seq.getTranslatedAA("AUG");
     => "M"
     seq.getTranslatedAA("GCU");
     => "A"
     seq.getTranslatedAA("CUU");
     => "L"
 ### Remove introns
 Take a sequence and an array of exonsRanges and removes them.

     seq.removeIntrons("ATGACCCTGAAGGTGAATGACAG", [[1, 8]]);
     => "TGACCCT"
     seq.removeIntrons("ATGACCCTGAAGGTGAATGACAG", [[2, 9], [12, 20]]);
     => "GACCCTGGTGAATGA"
 ### Transcribe sequence
 Takes a sequence string and returns the transcribed sequence (dna <-> rna).
 If an array of exons is given, the introns will be removed from the sequence.

     seq.transcribe("ATGACCCTGAAGGTGAA");
     => "AUGACCCUGAAGGUGAA"
     seq.transcribe("AUGACCCUGAAGGUGAA"); reverse
     => "ATGACCCTGAAGGTGAA"
 ### Translate sequence
 Takes a DNA or RNA sequence and translates it to protein
 If an array of exons is given, the introns will be removed from the sequence.

     seq.translate("ATGACCCTGAAGGTGAATGACAGGAAGCCCAAC"); dna
     => "MTLKVNDRKPN"
     seq.translate("AUGACCCUGAAGGUGAAUGACAGGAAGCCCAAC"); rna
     => "MTLKVNDRKPN"
     seq.translate("ATGACCCTGAAGGTGAATGACAGGAAGCC", [[3, 21]]);
     => "LKVND"
 ### Reverse exons
 Takes an array of exons and the length of the reference and returns inverted coordinates.

     seq.reverseExons([[2,8]], 20);
     => [ [ 12, 18 ] ]
     seq.reverseExons([[10,45], [65,105]], 180);
     => [ [ 135, 170 ], [ 75, 115 ] ]
 ### Find non-canonical splice sites
 Takes a sequence and exons ranges and returns an array of non canonical splice sites.

     seq.findNonCanonicalSplices("GGCGGCGGCGGTGAGGTGGACCTGCGCGAATACGTGGTCGCCCTGT", [[0, 10], [20, 30]]);
     => [ 20 ]
     seq.findNonCanonicalSplices("GGCGGCGGCGGTGAGGTGAGCCTGCGCGAATACGTGGTCGCCCTGT", [[0, 10], [20, 30]]);
     => []
 ### Check canonical translation start site
 Takes a sequence and returns boolean for canonical translation start site.

     seq.checkCanonicalTranslationStartSite("ATGACCCTGAAGGT");
     => true
     seq.checkCanonicalTranslationStartSite("AATGACCCTGAAGGT");
     => false
 ### Get reading frames
 Takes a sequence and returns an array with the six possible Reading Frames (+1, +2, +3, -1, -2, -3).

     seq.getReadingFrames("ATGACCCTGAAGGTGAATGACAGGAAGCCCAAC");
     => [ 'ATGACCCTGAAGGTGAATGACAGGAAGCCCAAC',
          'TGACCCTGAAGGTGAATGACAGGAAGCCCAAC',
          'GACCCTGAAGGTGAATGACAGGAAGCCCAAC',
          'GTTGGGCTTCCTGTCATTCACCTTCAGGGTCAT',
          'TTGGGCTTCCTGTCATTCACCTTCAGGGTCAT',
          'TGGGCTTCCTGTCATTCACCTTCAGGGTCAT' ]
 ### Get open reading frames
 Takes a Reading Frame sequence and returns an array of Open Reading Frames.

     seq.getOpenReadingFrames("ATGACCCTGAAGGTGAATGACAGGAAGCCCAAC");
     => [ 'ATGACCCTGAAGGTGAATGACAGGAAGCCCAAC' ]
     seq.getOpenReadingFrames("AUGACCCUGAAGGUGAAUGACAGGAAGCCCAAC");
     => [ 'AUGACCCUGAAGGUGAAUGACAGGAAGCCCAAC' ]
     seq.getOpenReadingFrames("ATGAGAAGCCCAACATGAGGACTGA");
     => [ 'ATGAGAAGCCCAACATGA', 'GGACTGA' ]
 ### Get all open reading frames
 Takes a sequence and returns all Open Reading Frames in the six Reading Frames.

     seq.getAllOpenReadingFrames("ATGACCCTGAAGGTGAATGACA");
     => [ [ 'ATGACCCTGAAGGTGAATGACA' ],
          [ 'TGA', 'CCCTGA', 'AGGTGA', 'ATGACA' ],
          [ 'GACCCTGAAGGTGAATGA', 'CA' ],
          [ 'TGTCATTCACCTTCAGGGTCAT' ],
          [ 'GTCATTCACCTTCAGGGTCAT' ],
          [ 'TCATTCACCTTCAGGGTCAT' ] ]

 ### Find longest open reading frame
 Takes a sequence and returns the longest ORF from all six reading frames and
 corresponding frame symbol (+1, +2, +3, -1, -2, -3). If a frame symbol is specified,
 only look for longest ORF on that frame.
 When sorting ORFs, if there's a tie, choose the one that starts with start codon Methionine.
 If there's still a tie, return one randomly.

     seq.findLongestOpenReadingFrame("ATGACCCTGAAGGTGAATGACA");
     => [ 'ATGACCCTGAAGGTGAATGACA', '+1' ]
     seq.findLongestOpenReadingFrame("ATGACCCTGAAGGTGAATGACA", "-1");
     => "TGTCATTCACCTTCAGGGTCAT"
