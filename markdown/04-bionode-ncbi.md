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
