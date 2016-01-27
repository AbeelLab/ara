// Written by Thies Gehrmann thiesgehrmann@gmail.com 2015-2016

/******************************************************************************
 * Determine Unique sequence identifiers based on on two inputs               *
 *                                                                            *
 ******************************************************************************/

package ara.uniquemarkers




object UniqueMarkers {

  /////////////////////////////////////////////////////////////////////////////

  def main(args: Array[String]) {

    if(args.length < 4) {
      usage("uniqueMarkers")
      System.exit(1);
    }

      // Read input
    val inputFastas  = args(0).split(',').toList
    val inputGenomes = args(1).split(',').toList
    val k            = args(2).toInt
    val out_prefix   = args(3)
    val redundant    = if (args.length > 4 && args(4) == "R"){ true } else { false }

      // Read the genes from different alleles
    val fastas  = inputFastas.map(Fasta.read);
      // Read genome FASTA files
    val genomes = inputGenomes.map(Fasta.read);

      // Determine the unique kmers for each group of alleles
    val orthologous_kmers = Utils.suckLessZip(fastas).map( o => KmerTools.orthologUniqueKmers(o, k))
      // Create the HashMaps for the genomes
    val genomeKmerCounts  = new genomeKmers(genomes.flatten, k);

      // Determine the genome-wide unique markers for each allele
    val uniqueMarkers = KmerTools.orthologGenomeUniqueKmers(orthologous_kmers, genomeKmerCounts)

      // Find nonredundant markers.
      // If requested, skip this step
    val outMarkers = if (redundant) { uniqueMarkers } else { KmerTools.nonRedundantKmers(uniqueMarkers, k) }

      // Output
    output(outMarkers, inputFastas.length, out_prefix)

    System.exit(0);

  }

  /////////////////////////////////////////////////////////////////////////////

  import java.io._

  def output(markers: List[(List[String], List[Set[kmerType]])], n_genomes: Int, prefix: String): Unit = {

    val outfd_fasta_all = new PrintWriter(new FileWriter(prefix + "_all.fasta", false))

    for(genome_i <- 0 until n_genomes) {

     val outfd = new PrintWriter(new FileWriter(prefix + '_' + genome_i.toString + ".tsv", false))
     val outfd_fasta = new PrintWriter(new FileWriter(prefix + '_' + genome_i.toString + ".fasta", false))

     println("--------Unique Markers for genome " + genome_i.toString + "--------")

      for(gene_i <- 0 until markers.length){
        val gene_markers = markers(gene_i)._2(genome_i).toList
        for(probe_i <- gene_markers.sortWith(_.index < _.index)) {
          println(markers(gene_i)._1(genome_i) + '\t' + probe_i.index.toString + '\t' + probe_i.seq)
          outfd.write(markers(gene_i)._1(genome_i) + '\t' + probe_i.index.toString + '\t' + probe_i.seq + '\n')
          outfd_fasta.write(">" + markers(gene_i)._1(genome_i) + '_' + probe_i.index.toString  + '\n');
          outfd_fasta_all.write(">" + markers(gene_i)._1(genome_i) + '_' + probe_i.index.toString  + '\n');
          for ( l <- probe_i.seq.grouped(80).toList) {
            outfd_fasta.write(l + '\n');
            outfd_fasta_all.write(l + '\n');
          }

        }
      }
      outfd.close()
      outfd_fasta.close()
    }
    outfd_fasta_all.close()

  }

  /////////////////////////////////////////////////////////////////////////////*/

  def usage(arg0: String): Unit = {
    println("Unique identifier detection")
    println("Usage: unique-markers  <sequence1>[,sequence2[,...]] <genome1>[,genome2[,...]] <k> <out_prefix> [R]");
    println("  sequences:  Comma seperated list of FASTA files containing gene DNA sequences")
    println("  genomes:    Comma seperated list of FASTA files containing genome DNA sequences")
    println("  k:          The size of kmer to use")
    println("  out_prefix: The prefix of the output files. Will produce *.tsv files for each input file")
    println("  R:          If this flag is present, the tool will produce redundant markers")
    println("")
    println("Example usage")
    println("  $>" + arg0 + " interestingGenes.fasta genomeSequence.fasta 21 interestingGeneMarkers")
    println("    For each gene listed in <interestingGenes.fasta>, the tool will find markers of length 21bp that are unique to that gene w.r.t. the sequences present in <genomeSequence.fasta>.")
    println("")
    println("  $>" + arg0 + " interestingGenes_org1.fasta,interestingGenes_org2.fasta, genome_org1.fasta,genome_org2.fasta 21 orthologMarkers")
    println("    For each gene pair of genes in interestingGenes_org1.fasta and interestingGenes_org2.fasta, the tool will find markers that distinguish the two pairs, and are unique w.r.t genome_org1.fasta and genome_org2.fasta.")
    println("")
    println("  $>" + arg0 + " interestingGenes.fasta genome_org1.fasta,genome_org2.fasta 21 interestingGeneMarkers R")
    println("    For each gene in interestingGenes.fasta, the tool will find markers that distinguish the gene from all the other sequences in genome_org1.fasta and genome_org2.fasta, but these markers will be overlapping.")

  }

  /////////////////////////////////////////////////////////////////////////////*/

}
