# Introduction

## What is [Bionode](http://bionode.io)?
It's "Modular and universal bioinformatics"! Modular because it follow the UNIX philosophy of small tools, each one just trying to do "one thing well". Universal because it can be used anywhere, on a laptop, a server or even a browser!

There is fundamentally three things that sets it apart from other [bio* libraries](http://www.open-bio.org/wiki/Projects#Main_projects):
* It can run in a browser since it's written in JavaScript, but also on a server since it's using [Node.js](http://nodejs.org);
* It takes advantage of [Node.js' asynchronous and event-based nature](https://nodejs.org/en/about), and [Streams](http://joshondesign.com/2014/06/25/nodestreamsareawesome) to pipe data or objects around in a scalable way;
* In addition to a [JavaScript API](https://github.com/bionode/bionode-ncbi#usage), each module also provides a [command line interface](https://github.com/bionode/bionode-ncbi#command-line-examples) that allows users to use and mix bionode in their existing tool chain of favorite programming language, without having to learn or write JavaScript.

Each module has it's own [GitHub repository](http://github.com/bionode) with testing and documentation, but there's a [main module](http://github.com/bionode/bionode) that simply depends on the other modules, so that they can all be installed in one step.

Modules exchange data with each other using the standards [NDJSON](http://ndjson.org) and [Protocol Buffers](https://developers.google.com/protocol-buffers).

## How did it started?
[Bruno Vieira](http://bmpvieira.com) started it in 2014 as a side project of his PhD in Bioinformatics and Population Genomics at [WurmLab](http://wurmlab.github.io) because he needed to write code for several purposes, and some of it would eventually end up having to be rewritten in JavaScript for some [biological web projects](http://wurmlab.github.io/tools/).
Since he had a background in Node.js, he decided to just write everything in Node.js from the start, and also benefit of some of Node.js' features that seemed very appropriate for orchestrating bioinformatics pipelines. The philosophy of small modules that is present in the Node.js community also allowed for each piece of code to be quickly published as a module with a few extras (continous integration, auto-generated documentation, badges) rather than getting lost as scripts somewhere waiting to be integrated in a bigger monolithic project.

[![badges](/static/img/badges.png)](http://github.com/bionode/bionode-ncbi)

## How big is the community?
Bionode has a small community around it, with a few active developers that already contributed new modules.

![team](/static/img/team.png)

Bionode's best practices have also inspired [other projects](https://github.com/hydronode/hydronode).

## What kind of tools are you developing?
The overall theme is bioinformatics, besides that, any tools is welcome. The idea is that if someone is writing JavaScript for biology, or agrees with some of our [principles](http://github.com/bionode/bionode-template#principles), then that code is welcome in Bionode.
Currently there's modules to retrieve biological data from web resources (e.g., [bionode-ncbi](http://github.com/bionode/bionode-ncbi)), deal with specific bioinformatic data formats (e.g., [bionode-fasta](http://github.com/bionode/bionode-fasta)), and handle sequences ([bionode-seq](http://github.com/bionode/bionode-seq)). There's also wrappers around existing bioinformatics tools (e.g., [bionode-sam](http://github.com/bionode/bionode-sam)), which obviously only run server side (since those tools are not written in JavaScript) to make it easier to build [pipelines in Node.js](https://github.com/bionode/bionode-example-dat-gasket#bionode-example-with-dat-and-gasket). However there is some discussion if [wrappers are in the scope](http://github.com/bionode/bionode/issues/29) of Bionode.
Once we have more basic tools, there might be a move towards genetic algorithms and machine learning.

## Current state
You can see a global list of modules on [GitHub](https://github.com/bionode/bionode#list-of-modules). Many are still in early development, so check the READMEs and the status badges. On each README you have a [link to the documentation](https://github.com/bionode/bionode-ncbi#usage) of each module, but you can also check a concatenated version of all the documentation at http://doc.bionode.io.
You can also have an overview of all the issues at [waffle.io](http://waffle.io/bionode/bionode). Theres plenty of opportunity for collaboration, and contributors will be coauthors on a paper, so check the ["How to contribute"](/guide/11-how-to-contribute.html) section.  

## Get involved
You can also check our [Twitter account](http://twitter.com/bionode) for updates or come ask questions in the [Gitter chat room](http://gitter.im/bionode/bionode). You can also just open an [issue on GitHub](http://github.com/bionode/bionode/issues) with your question or idea. Have a look around and come chat to us, there's plenty of ways to get involved, from either the biology side or the computational side.

## Collaborations
### Dat
Bionode collaborates with [Dat](http://dat-data.com), a "Git for data". The Dat developers have contributed and helped bionode, and both projects follow similar standards. Bionode influenced some of the development decisions in Dat, towards the goal of making Dat the definitive format to share scientific data.

![sanger dat/bionode talk](/static/img/sanger.jpg)

### BioJS
[BioJS](http://biojs.net) is a set of libraries for vizualising biological data on the web. It also provides as [registry of biological tools](http://biojs.io) on the web. Bionode and BioJS collaborate to improve compatibility and integration between both projects, for example, by using the [same package manager (npm)](http://github.com/bionode/bionode/issues/9). We also submit [Google Summer of Code projects](http://biojs.net/gsoc/2015/) together.

Since both projects are JavaScript and have the word "bio" in the name, people ask frequently what is the difference between them? There are several, but they can be summarized by the fact that Bionode is focusing and "data processing/pipelines" (more backend, even when in a browser) while BioJS e focusing on "visualization". This means that even if the language is the same (which is great for integration) the architecture of each project is very different. See the Venn diagram below with some of the technologies used by each project:

![venn](/static/img/venn.png)

## Node.js? Really??
JavaScript is currently the only real "write once, run anywhere" language, and for most cases it's fast enough. Nowadays you can even run very [demanding 3D apps](https://blog.mozilla.org/blog/2014/03/12/mozilla-and-epic-preview-unreal-engine-4-running-in-firefox/) compiled to JavaScript in you're browser. If you want, you can also run [C code in Node.js for speed](http://www.benfarrell.com/2013/01/03/c-and-node-js-an-unholy-combination-but-oh-so-right/).
The asynchronous and evented nature of Node.js also makes if very interesting for big data pipelines.

Plus it's the fastest growing open source community!

![modules](/static/img/modules.png)

Sure, it still doesn't have has many bio tools as Python or statistics as R, but that's an opportunity to contribute and get involved in some of the other cool things the Node.js community is making.

## Talks and workshops

* [Dat workshop with Bionode chapter (8) at MozFest 2014 London, UK](http://try-dat.com)
* [Dat and Bionode at Sanger, UK](https://www.youtube.com/watch?v=Ef17lkx7s0U)
* [Open Research Cambridge meetup, UK](http://www.eventbrite.co.uk/e/building-collaborative-workflows-for-scientific-data-tickets-14527561327)
* [Biohackathon in Nagasaki, Japan](https://www.youtube.com/watch?v=9MoI1IFRdvc&index=15&list=PL0uaKHgcG00bSajcVd8qIQ__Mss8xkPoH)
* [Biocoders meetup in Cambridge, UK](http://www.meetup.com/biocoders/events/225520856/)
* MozFest 2015 London, UK
