# Getting started
## Online (this workshop)
If you're planning to try Bionode out using [try.bionode.io](http://try.bionode.io), you're already all set!

## Install it on your machine (alternative)
Otherwise, you want to run it locally on you machine, the fist step is to get Node.js. There are several ways to do it, but we recommend the following:

```bash
# OS X (using Homebrew)
brew install n # Get a Node.js version manager, very useful.

# Ubuntu
sudo apt-get install npm
npm install n -g

# Windows
## Go to http://nodejs.org, and follow instructions

# Install version 0.10.40 of Node.js and load it in your $PATH
n 0.10.40

# Get Bionode and friends
## Installs bionode 'globally', i.e., as a Command Line
npm install bionode -g
## You can also install a specific module instead of all
npm install bionode-ncbi -g
## Install some other useful tools
npm install dat json -g
```

## Available modules
After you're setup you can have a quick look at the available [modules on GitHub](https://github.com/bionode/bionode#list-of-modules) and jump to the section about that module, or keep reading

## How things work in general

### Command Line Interface
Check the documentation and status of each module, but in general you can use the command line interface like this:

```bash
bionode ncbi urls assembly Acromyrmex | json -ga genomic.fna
```

That command queries the NCBI database and retrieves URLs of the genome assembly for the ant species Acromyrmex. This will return a JSON object that is then piped to the `json` command so that we can retrieve only the property `genomic.fna` (the url of the file with DNA sequences in fna/fasta format) and filter out the other properties.

### JavaScript API
Now the same could be done using the JavaScript API, but first you need to create a folder for your project and then for each module your are going to `require` in your code, you need to do `npm install module_name` (without the `-g` flag) to install a copy of that module locally in your project folder. You only use the `-g` flag when you want to install a module as a command line tool.

```bash
npm install bionode-ncbi
```

```javascript
var bio = require('bionode')
```

 The you can generally use bionode modules in 3 different ways:

#### The Callback pattern
A [callback](https://docs.nodejitsu.com/articles/getting-started/control-flow/what-are-callbacks) simply means, you ask for something and once you get all of it you process it

```javascript
// Query NCBI
bio.ncbi.urls('assembly', 'Acromyrmex', function(urls) {
  # Got all the urls as an array, print just first genome
  console.log(urls[0].genomic.fna)
})
```

#### The Event pattern
Callbacks are fine for most cases, but if you're getting too much data your code will run out memory and crash. A solutions is to use [Events](https://nodesource.com/blog/understanding-the-nodejs-event-loop) to do something as you get one object or chunks of data.

```javascript
bio.ncbi.urls('assembly', 'Acromyrmex').on('data', printGenomeURL)
function printGenomeURL(url) {
  console.log(url.genomic.fna)
}
```

#### The Pipe pattern
[Node.js Streams](https://github.com/substack/stream-handbook) are based on Events and allow you to get rid of a lot of boilerplate code by chaining functions together.

```javascript
var tool = require('tool-stream')
bio.ncbi.urls('assembly', 'Acromyrmex')
.pipe(tool.extractProperty('genomic.fna'))
.pipe(process.stdout)
```

### How is it done in other libraries?

Here's an example of how you would do the same in [BioPython](http://biopython.org) (other libs are similar):

```python
# URL for the Acromyrmex assembly?
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000188075.1_Si_gnG
import xml.etree.ElementTree as ET
from Bio import Entrez
Entrez.email = "mail@bmpvieira.com"
esearch_handle = Entrez.esearch(db="assembly", term="Acromyrmex")
esearch_record = Entrez.read(esearch_handle)
for id in esearch_record['IdList']:
  esummary_handle = Entrez.esummary(db="assembly", id=id)
  esummary_record = Entrez.read(esummary_handle)
  documentSummarySet = esummary_record['DocumentSummarySet']
  document = documentSummarySet['DocumentSummary'][0]
  metadata_XML = document['Meta'].encode('utf-8')
  metadata = ET.fromstring('' + metadata_XML + '')
  for entry in Metadata[1]:
    print entry.text
```

### More Node.js tips

If you git clone a Node.js folder, to install it's dependencies you just cd into it and type `npm install`.
If you want to install that module that you just git cloned as a command tool, you cd into the folder and do `npm link` (useful for development).
