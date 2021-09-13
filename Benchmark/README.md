# Correlation Coefficent Benchmark for Genetic Recombination
## Instructions
`./compCC.sh <resolution> <genetic map>`

Example: `./compCC.sh 10000 output.map`

## Genetic Map Format
The format of the genetic map for every line must be:
`<genomic position> <recombination rate> <genetic mapping>`

The first line of the genetic map must be a header.

## Output Format
The output will be print on System output and will consist of 4 space seperated integers:
`<2M> <5M> <midregion> <overall>`

