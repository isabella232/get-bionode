# Write Once Run Anywhere (WIP section)

Node.js is JavaScript, and thus can also run on a browser. However the Node.js API won't be available and things like filesystem access don't make sense in the browser, so the workaround is to use the project [browserify](http://browserify.org) that makes a Node.js module work in the browser by implementing the Node.js API with browser functions. For example, filesystem works by storing files in the browser leveldb database.

The module [bionode-fasta](https://github.com/bionode/bionode-fasta) is a good example of WORA code. You can look at its test folder to see what's going to happen and then type ```npm run test-browser``` as shown below to see it run those test in a browser.


```bash
git clone https://github.com/bionode/bionode-fasta
cd bionode-fasta
npm install
npm run test-browser
```

The result should be this:

![client-side](/static/img/client-side.png)
