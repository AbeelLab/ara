package macaw2

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object QCsnpFasta {

  case class Config(
    val fastaFile: File = null,
    val output: File = null,
    val function: String = null,
    val threshold: Double = 0.05
  )

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar qc-snp-fasta") {
      opt[File]('i', "infile") required () action { (x, c) => c.copy(fastaFile = x) } text ("Fasta file with SNP sequences.")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output file.")
      arg[String]("<samples | N-content>") required () validate { x => if (x == "samples" || x == "N-content") success else failure("Required argument <samples | N-content>") } action { (x, c) => c.copy(function = x) } text ("Choose function.")
      opt[Double]('t', "threshold") validate { x => if (x.isInstanceOf[Double]) success else failure("Wrong value")} action { (x, c) => c.copy(threshold = x) } text ("Default threshold for N-content = 0.05")
    }

    def getSamples(src: Source): List[String] = {
      src.getLines.filter { _.startsWith(">") }.toList
    }

    def getLength(src: Source): (String, Int) = {
      val s = src.getLines.take(1).mkString
      val seq = src.getLines.takeWhile { x => !x.startsWith(">") }.mkString
      (s, seq.size)
    }

    def getNcontent(src: Source, t: Double, pw: PrintWriter): Unit = {
      pw.println("#Sample\tN-count\tN-percentage\t>threshold\tTotal Length")
      var it = src.getLines
      var highNsamples = 0
      while (!it.isEmpty) {
        var s = it.take(1).mkString
        var (seq, res) = it.span { x => !x.startsWith(">") }
        var seqStr = seq.mkString
        var Ncount = seqStr.count { x => x == 'N' }
        var Npercentage = Ncount.toFloat / seqStr.size
        pw.println(s + "\t" + Ncount + "\t" + Npercentage + "\t" + (if (Npercentage > 0.05) 1 else 0) + "\t" + seqStr.size)
        if (Npercentage > 0.05) {
          highNsamples += 1
        }
        it = res
      }
      pw.println("#")
      pw.println("# " + highNsamples + " samples with high N content (>0.05%).")
    }

    parser.parse(args, Config()) map { config =>
      val pw = new PrintWriter(config.output)
      val src = Source.fromFile(config.fastaFile)
      val t = config.threshold
      
      config.function match {
        case "N-content" => getNcontent(src, t, pw)
        case "samples" => getSamples(src)
      }
      pw.close
    }
    
    
  }
}