package ara.dr

import scala.io.Source
import java.io.File
import java.io.PrintWriter
import ara.Gene._
import ara.dr.DRsnp._
import ara.Gene

/**
 * Create minimized version of reference genome H37Rv that consists of genes and their flanking regions containing drug resistances from the list.
 *
 */
object DRMap {

  case class Config(
    val drList: File = null,
    val fasta: File = null,
    val output: String = null,
    val gff: File = null
  )
    
    
  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar ref-DR-region") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(drList = x) } text ("List of TB drug resistances.")
      opt[File]("ref") required () action { (x, c) => c.copy(fasta = x) } text ("Reference fasta file.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("name of output file.")
      opt[File]("gff") required () action { (x, c) => c.copy(gff = x) } text ("gff-file corresponding to reference fasta-file.")
    }

    parser.parse(args, Config()) map { config =>

      /** Configure files */
      val ref = Source.fromFile(config.fasta).getLines
      val refName = ref.take(1).mkString
      val refGenome = ref.filterNot(_.startsWith(">")).mkString
      val pw = new PrintWriter(new File(config.output))

      /** Read SNPs from DR list */
      val drList = Source.fromFile(config.drList).getLines.filterNot(_.startsWith("#")).filter(_.isSNP).map { line =>
        line match {
          case DRsnp(d, l, lt, cp, r, gp, a, cn, aac) => new DRsnp(d, l, lt, cp, r, gp, a, cn, aac)
        }
      }.toList

      /** Map locus info */
      val allGenes = drList.map { g =>
        val locus = if (g.locus.endsWith("promoter")) g.locus.split("-")(0) else g.locus
        (g.locusTag, locus)
      }.distinct.sorted
      //allGenes.foreach(println)
      //println("\nallGenes.size: " + allGenes.size + "\n")

      /** Read gff file and filter only genes that occur in DR list. */
      val gff = Source.fromFile(config.gff).getLines.filterNot { line =>
        val feature = line.split("\t")(2)
        feature == "CDS" || feature == "source"
      }.map {
        _ match {
          case Gene(tag, feature, start, end, dir) => (tag -> new Gene(tag, feature, start, end, dir))
        }
      }.toMap
      //gff.toList.sortBy(_._1).foreach(println)
      //println("\ngff.size: " + gff.size)

      /** Get Gene info for each locus, and group genes together if they are within 1000 bp of each other. */
      def ranges(ls: List[(Gene, String)], range: Int): List[List[(Gene, String)]] = ls match {
        case head :: tail => {
          var (b, a) = ls.span(p => head._1.end + range > p._1.start - range)
          if (b.size > 1) {
            while (!a.isEmpty && b.last._1.end + range > a.head._1.start - range) {
              b = b :+ a.head
              a = a.drop(1)
            } 
          }
          b :: ranges(a, range)
        }
        case Nil => Nil
      }

      
      val genes = allGenes.map(g => (gff(g._1), g._2)).sortBy(_._1.start)
      val allRanges = ranges(genes, 1000)
      allRanges.foreach(println)
      println("\nallRanges.size: " + allRanges.size)
      println("\n" + allRanges.flatMap(g => g).size + " genes")
   
      /** Print multi-fasta file, each gene and 1000 bp flanking regions. */
     allRanges.foreach { g =>
        val start = g.head._1.start - 1000
        val end = g.last._1.end + 1000
        val genes = g.map(_._2).mkString("/")
        pw.println(refName + "_" + genes + "_" + start + "-" + end)
        refGenome.substring(start - 1, end).grouped(80).foreach(pw.println)
        pw.println
      }

      pw.close
    }

  }
}