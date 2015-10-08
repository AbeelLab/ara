package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import be.abeel.bioinformatics.FastaIterator

object QCsnpFasta {

  case class Config(
    val fastaFile: File = null,
    val output: File = null,
    val function: String = null,
    val vcfPathFile: File = null,
    val dirIdx: Int = 0
  )

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar qc-snp-fasta") {
      opt[File]('i', "infile") required () action { (x, c) => c.copy(fastaFile = x) } text ("Fasta file with SNP sequences.")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output file.")
      arg[String]("<samples | N-content>") required () validate { x => if (x == "samples" || x == "N-content") success else failure("Required argument <samples | N-content>") } action { (x, c) => c.copy(function = x) } text ("Choose function.")
      opt[File]("path-file") required () action { (x, c) => c.copy(vcfPathFile = x) } text ("Remove samples with high N-content from VCF path file.")
      opt[Int]("index") required () action { (x, c) => c.copy(dirIdx = x) } text ("Directory index of sample name in VCF path file.")
    }

    def getSamples(src: Source): List[String] = {
      src.getLines.filter { _.startsWith(">") }.toList
    }

    def getLength(src: Source): (String, Int) = {
      val s = src.getLines.take(1).mkString
      val seq = src.getLines.takeWhile { x => !x.startsWith(">") }.mkString
      (s, seq.size)
    }

    def getNcontent(src: File, vcfPaths: Map[String, String], pw1: PrintWriter, pw2: PrintWriter): Unit = {
      pw1.println("#Sample\tN-count\tN-percentage\t>threshold\tTotal Length")
      var highNsamples = 0
      val fit = new FastaIterator(src)
      while (fit.hasNext()) {
        val item = fit.next()
        val header = item.getDescription.drop(1)
        val seq = item.getSequence
        val Ncount = seq.count { _ == 'N' }
        var Npercentage = Ncount.toFloat / seq.size
        pw1.println(header + "\t" + Ncount + "\t" + Npercentage + "\t" + (if (Npercentage > 0.05) 1 else 0) + "\t" + seq.size)
        if (Npercentage > 0.05) {
          highNsamples += 1
        } 
        else if (header != "reference") {
          pw2.println(vcfPaths(header))
        }
      }
      pw1.println("#")
      pw1.print("# " + highNsamples + " samples with high N content (>0.05%).")
    }

    /** Elapsed time function */
    def time[R](block: => R): R = {
      val t0 = System.currentTimeMillis()
      val result = block // call-by-name
      val t1 = System.currentTimeMillis()
      println("Elapsed time: " + (t1 - t0) + "ms")
      result
    }

    parser.parse(args, Config()) map { config =>
      val pw1 = new PrintWriter(config.output)
      val src = Source.fromFile(config.fastaFile)
      val pathFile = config.vcfPathFile
      val pw2 = new PrintWriter(pathFile.getName + "_QC")
      val vcfPaths = Source.fromFile(pathFile).getLines.map(p => (p.split("/")(config.dirIdx), p)).toMap
      time {
        config.function match {
          case "N-content" => getNcontent(config.fastaFile, vcfPaths, pw1, pw2)
          case "samples" => getSamples(src)
        }
        pw1.close
        pw2.close
      }

    }

  }
}