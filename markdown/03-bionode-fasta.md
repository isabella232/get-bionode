# bionode-fasta
> Streamable FASTA parser.
>
> doi: [?](?)
> author: [Bruno Vieira](http://bmpvieira.com)
> email: <mail@bmpvieira.com>
> license: [MIT](https://raw.githubusercontent.com/bionode/bionode-fasta/master/LICENSE)

---

## Usage
This module can be used in Node.js as described further below, or as a command line tool.
Examples:

    $ npm install -g bionode-fasta

    # bionode-fasta [options] [input file] [output file]
    $ bionode-fasta input.fasta.gz output.json

    # You can also use fasta files compressed with gzip
    # If no output is provided, the result will be printed to stdout
    # Options: -p, --path: Includes the path of the original file as a property of the output objects
## Fasta
Returns a Writable Stream that parses a FASTA content Buffer into a JSON Buffer

    var fasta = require('bionode-fasta')

    fs.createReadStream('./input.fasta')
    .pipe(fasta())
    .pipe(process.stdout)

    => { "id": "contig1",
         "seq": "AGTCATGACTGACGTACGCATG" }
    => { "id": "contig2",
         "seq": "ATGTACGTACTGCATGC" }
    => [...]

Can also parse content from filenames Strings streamed to it

    fs.createReadStream('./fasta-list.txt')
    .pipe(split())
    .pipe(fasta({filenameMode: true}))
    .pipe(process.stdout)

When filenames are Streamed like in the previous example, or passed directly
to the parser Stream, they can be added to the output Objects

    fasta({includePath: true}, './input.fasta')
    .pipe(process.stdout)

    => { "id": "contig1",
         "seq": "AGTCATGACTGACGTACGCATG" }
         "path": "./input.fasta" }

The output from the parser can also be available as Objects instead of Buffers

    fasta({objectMode: true}, './input.fasta')
    .on('data', console.log)

Shortcut version of previous example

    fasta.obj('./input.fasta').on('data', console.log)

Callback style can also be used, however they might not be the best for large files

    fasta.obj('./input.fasta', function(data) {
      console.log(data)
    })