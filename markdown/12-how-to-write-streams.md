# How to write Node.js Streams (WIP section)
Update (#MozFest): Check the [mississippi module](https://github.com/maxogden/mississippi).

```javascript
var through = require('through2')
var stream = through2.obj(transform)
function transform (obj, enc, next) {
  // do things, example:
  obj.name = obj.name.toUpperCase()
  // Push downstream
  this.push(obj)
  // Callback to fetch next object
  next()
}

var through = require('through2')
var stream = through2.obj(transform)
function transform (obj, enc, next) {
  // do things, example:
  var self = this
  requestSomethingFromDB(obj.name, function(data) {
    obj.data = data
    self.push(obj)
    next()
  })
}
```

```bash
mkdir project
cd project
npm install bionode-ncbi through2
```

```javascript
var ncbi = require('bionode-ncbi')
var through = require('through2')
var json = require('ndjson')

var myStream = through.obj(transform)
function transform (obj, enc, next) {
  var result = {
    specie: obj.organism,
    organisazation: obj.meta['submitter-organization']
  }
  this.push(result)
  next()
}

ncbi.search('assembly', 'spiders')
.pipe(myStream)
.pipe(json.stringify())
.pipe(process.stdout)
```


```javascript
var counter = 0
  myStream
  .on('data', function (data) {
    counter++
  })
  .on('end', function () {
    console.log('Processed ' + counter)
  })


  var counter = 0

  var count = function (data) {
    counter++
  }

  var log = function () {
    console.log('Processed ' + counter)
  }

  myStream.on('data', count).on('end', log)
```
