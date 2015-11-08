# Bionode-ncbi (WIP section)

```bash
bionode-ncbi search genome spiders
bionode-ncbi search genome spiders | wc
bionode-ncbi search genome spiders | head -n 1 | json
bionode-ncbi search genome spiders | json -ga organism_name

bionode-ncbi search genome spiders | \
  json -ga uid | \
  bionode-ncbi link genome pubmed - | \
  json -ga destUID | \
  bionode-ncbi search pubmed - | \
  json -ga title

bionode-ncbi download assembly Guillardia theta | \
  json -ga -c 'this.status === "completed"' | \
  json -ga path | \
  bionode-fasta -f | \
  json -ga -c 'this.seq.length > 10000' | \
  bionode-fasta --write > gtheta-big-scaffolds.fasta
```

# bionode-ncbi
> Node.js module for working with the NCBI API (aka e-utils) using Streams.
>
> doi: [10.5281/zenodo.10610](http://dx.doi.org/10.5281/zenodo.10610)
> author: [Bruno Vieira](http://bmpvieira.com)
> email: <mail@bmpvieira.com>
> license: [MIT](https://raw.githubusercontent.com/bionode/bionode-ncbi/master/LICENSE)
//
---

//
## Usage
This module can be used in Node.js as described further below, or as a command line tool.
Examples:
//
    $ npm install -g bionode-ncbi
//
    # bionode-ncbi [command] [arguments] --limit (-l) --throughput (-t)
    $ bionode-ncbi search taxonomy solenopsis
    $ bionode-ncbi search sra human --limit 500 # only return 500 items
    $ bionode-ncbi search sra human --throughput 250 # fetch 250 items per API request
    $ bionode-ncbi download assembly solenopsis invicta
    $ bionode-ncbi urls sra solenopsis invicta
    $ bionode-ncbi link assembly bioproject 244018
    $ bionode-ncbi search gds solenopsis | dat import --json

## Search
Takes a NCBI database string and a optional search term and returns a stream of objects found:

    ncbi.search('sra', 'solenopsis').on('data', console.log)
    => { uid: '280116',
         expxml: {"Summary":{"Title":"Single Solenopsis invicta male","Platform":{"_":"ILLUMINA", [...],
         runs: {"Run":[{"acc":"SRR620577","total_spots":"23699662","total_bases":"4787331724", [...],
         extlinks: '    ',
         createdate: '2013/02/07',
         updatedate: '2012/11/28' }
    => { uid: '280243',
         expxml: {"Summary":{"Title":"Illumina small-insert paired end","Platform":{"_":"ILLUMINA", [...],
         runs: {"Run":[{"acc":"SRR621118","total_spots":"343209818","total_bases":"34320981800", [...],
         extlinks: '    ',
         createdate: '2013/02/07,
         updatedate: '2012/11/28' }
    => [...]

Arguments can be passed as an object instead:

    ncbi.search({ db: 'sra', term: 'solenopsis' })
    .on('data', console.log)

Advanced options can be passed using the previous syntax:

    var options = {
       db: 'assembly', // database to search
       term: 'human',  // optional term for search
       limit: 500,     // optional limit of NCBI results
       throughput: 100 // optional number of items per request
     }

The search term can also be passed with write:

     var search = ncbi.search('sra').on('data', console.log)
     search.write('solenopsis')

Or piped, for example, from a file:

     var split = require('split')

     fs.createReadStream('searchTerms.txt')
     .pipe(split())
     .pipe(search)


## Link
 Takes a string for source NCBI database and another for destination db and returns
 a objects stream with unique IDs linked to the passed source db unique ID.

     ncbi.link('taxonomy', 'sra', 443821)
     => { "srcDB":"taxonomy",
          "destDB":"sra",
          "srcUID":"443821",

          "destUID":"677548" }
     => { "srcDB":"taxonomy",
          "destDB":"sra",
          "srcUID":"443821",
          "destUID":"677547" }
     => [...]

 Also works with write and pipe, like **Search**.

## Property link (Plink)
 Similar to Link but taked the srcID from a property of the Streamed object
 and attached the result to a property with the name of the destination DB.

     ncbi.search('genome', 'arthropoda')
     .pipe(ncbi.expand('tax'))
     .pipe(ncbi.plink('tax', 'sra')

## Download
 Takes a NCBI database string and a optional search term and downloads the datasets/sequence files.
 ** Currently only supports sra and assembly databases. **
 Also accepts the keyword gff for annotations.
 Returns a stream that emits download progress and ends with download path
 The name of the folder where the file is saved corresponds to the UID from NCBI.

     ncbi.download('assembly', 'solenopsis invicta')
     .on('data', console.log)
     .on('end', function(path) {
       console.log('File saved at ' + path)
     }
     => Downloading 244018/unplaced.scaf.fa.gz 0.94 % of 106 MB at 0.48 MB/s
     => Downloading 244018/unplaced.scaf.fa.gz 100.00 % of 106 MB at 0.49 MB/s"
     => File saved at 244018/unplaced.scaf.fa.gz

## URLs
Takes a NCBI database string and a optional search term and returns as stream of dataset/sequence files URLs.
** Currently only supports sra and assembly databases. **
Also accepts the keyword gff for annotations.
The value of the uid property corresponds to the UID from NCBI.

     ncbi.urls('assembly', 'solenopsis invicta')
     .on('data', console.log)
     => {"url":"http://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/invertebrates/Solenopsis_invicta/Si_gnG/Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fa.gz",
         "uid":"244018/"}

## Expand
Takes a property (e.g., biosample) and optional destination property
(e.g., sample) and looks for a field named property+id (biosampleid)
in the Streamed object. Then it will do a ncbi.search for that id and save
the result under Streamed object.property.

    ncbi.search('genome', 'arthropoda').pipe(ncbi.expand('assembly'))

## Fetch
Allows retrieval of records from NCBI databases. Takes the database name, and a search term,
and returns the records from the database that match the search term. There are optional
advanced parameters that allow you to define how many records to retrieve and extra options
for genes. These parameters should be passed as an object.

i.e it can return a subset of a genetic sequence of a requested species

     ncbi.fetch('sra', 'solenopsis_invicta')
     => {"EXPERIMENT_PACKAGE_SET":
           {"EXPERIMENT_PACKAGE":
             [{"EXPERIMENT":
               [{"$":{"xmlns":"","alias":"Me","accession":"SRX757228,
               ...

With advanced parameters for sequence databases(all are optional):

     var opts = {
       db: 'nucest',
       term: 'guillardia_theta',
       strand: 1,
       complexity: 4
     }
     ncbi.fetch(opts)
     => { id: 'gi|557436392|gb|HE992975.1|HE992975 HE992975 Guillardia theta CCMP 327 Guillardia theta cDNA clone sg-p_014_h06, mRNA sequence',
         seq: 'GAAGGCGATTCCAATGGTGCGAGCGAGGCAGCGAACAGACGCAGCGGGGAGAG...
        }
     => [...]
For some databases there are multiple return types. A default one will be chosen
automatically, however it is possible to specify this via the rettype option e.g:

The NCBI website provides a list of databasese supported by efetch here:
http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly