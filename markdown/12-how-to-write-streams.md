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


 mkdir project
  cd project
  npm install bionode-ncbi through2



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
