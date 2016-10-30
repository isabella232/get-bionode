# Getting started

## TL; DR
Try this on the terminal on the right

```bash
bionode ncbi search nucleotide cancer | head -n 1 | json

bionode ncbi search nucleotide cancer | head -n 1 | json uid | bionode ncbi fetch nuccore - | json
```

## Prerequisites
We assume that you have some familiarity with a Command Line Interface (e.g., BASH).
If that is not the case, we recommend doing the [command_line_bootcamp](http://rik.smith-unna.com/command_line_bootcamp).
At a minimum, you need to know how to use the commands ```ls```, ```cd```, ```mkdir``` and ```touch```.

Knowledge of JavaScript and Node.JS is **not** required but can be very helpful for some sections. A good resource is [NodeSchool](http://nodeschool.io), and we recommend the sections ```javascripting```, ```learnyounode```, ```how-to-npm```, ```stream-adventure```, ```async-you```,  ```browserify-adventure```, and ```unctional-javascript-workshop``` (in this order). If you want a good and free JavaScript beginners book, check out [JavaScript for Cats](http://jsforcats.com).

If you get interested in Bioinformatics and want to learn more, there are plenty of resources and [MOOCs](https://en.wikipedia.org/wiki/Massive_open_online_course) out there. However, [Bioinformatics Data Skills](http://shop.oreilly.com/product/0636920030157.do) is a good beginners book.

## Do it online (this workshop)
You can test drive Bionode online without installing anything using the [try.bionode.io](http://try.bionode.io) website.

## Install it on your machine (alternative)
If you want to run it locally on your machine, the fist step is to get Node.JS. There are several ways to do it, but we recommend the following:

### Mac OS
Install Homebrew by copy pasting the following command in your terminal.

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Then install a Node.JS version manager
```bash
brew install n
```

Then install the latest stable version of Node or a specific one.

```bash
n stable
# Or
n 7.0.0
```

### Ubuntu
Run the following commands
```bash
# Install the Node Package Manager
sudo apt-get install npm
# Install a Node version manager
npm install n -g
# Install Node
n stable
## Or for a specific version
n 7.0.0
```
### Windows
Go to http://nodejs.org, and follow instructions.

## Install Bionode and other useful tools
Bionode provides a meta-module named ```bionode``` that can install all the other modules as dependencies. If you only need a specific module, you just install that one, e.g., ```bionode-ncbi```. Tip: In this tutorial, in the interest of  speed, saving computational resources, and avoid issues with some versions of Node.JS, we use ```--production``` to skip installing development dependencies.

### Installs bionode 'globally', i.e., as a Command Line tool (using ```-g```).
```bash
npm install bionode -g --production
```
### You can also install a specific module instead of all
```bash
npm install bionode-ncbi -g --production
```
### Install some other useful tools
```bash
npm install json tool-stream -g --production
```

## Available modules
After you're setup you can have a quick look at the available [modules on GitHub](https://github.com/bionode/bionode#list-of-modules) and jump to the section about that module, or keep reading

## How things work in general

### Command Line Interface
Check the documentation and status for each module in the README.md file on their GitHub page (e.g., [bionode-ncbi] (https://github.com/bionode/bionode-ncbi)), but in general you can use the command line interface like this:

```bash
bionode ncbi urls assembly Acromyrmex | json -ga genomic.fna
```

That command queries the NCBI database and retrieves URLs of the genome assembly for the ant species Acromyrmex. This will return a JSON object that is then piped to the `json` command so that we can retrieve only the property `genomic.fna` (the url of the file with DNA sequences in fna/fasta format) and filter out the other properties.

### JavaScript API
Now the same could be done using the JavaScript API, but first you need to create a folder for your project and then for each module your are going to `require` in your code, you need to do `npm install module_name` (without the `-g` flag) to install a copy of that module locally in your project folder. You only use the `-g` flag when you want to install a module as a command line tool.

```bash
#!/bin/bash
npm install bionode-ncbi
```

```javascript
#!/usr/bin/env node
var bio = require('bionode')
```

## Bionode code patterns
You can generally use bionode modules in 3 different ways:

### The Callback pattern
A [callback](https://docs.nodejitsu.com/articles/getting-started/control-flow/what-are-callbacks) simply means, you ask for something and once you get all of it you process it

```javascript
// Query NCBI
bio.ncbi.urls('assembly', 'Acromyrmex', function(urls) {
  # Got all the urls as an array, print just first genome
  console.log(urls[0].genomic.fna)
})
```

### The Event pattern
Callbacks are fine for most cases, but if you're getting too much data your code will run out memory and crash. A solutions is to use [Events](https://nodesource.com/blog/understanding-the-nodejs-event-loop) to do something as you get one object or chunks of data.

```javascript
bio.ncbi.urls('assembly', 'Acromyrmex').on('data', printGenomeURL)
function printGenomeURL(url) {
  console.log(url.genomic.fna)
}
```

### The Pipe pattern
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

## More Node.js tips

If you git clone a Node.js folder, to install its dependencies you can just cd into it and type `npm install`.
If you want to install that module that you just git cloned as a command line tool, you cd into the folder and do `npm link` (useful for development).
