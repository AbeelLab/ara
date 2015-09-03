package macaw2

import scala.io.Source
import com.sun.org.apache.xalan.internal.xsltc.compiler.Sort
import java.io.PrintWriter
import java.io.File
import macaw2.Mutation._

/**
 * This tool reads VCFs from the path file, and returns a SNP phy-file.
 */
object Vcf2snpPhylip {
  val usage = "java -jar ara.jar vcf2snp-phylip [pathfile] [output.phy]"

  def main(args: Array[String]) {

    /** Elapsed time function */
    def time[R](block: => R): R = {
      val t0 = System.currentTimeMillis()
      val result = block // call-by-name
      val t1 = System.currentTimeMillis()
      println("Elapsed time: " + (t1 - t0) + "ms")
      result
    }


    /** Get all SNP positions */
    def getPositions(fList: List[File]): Map[Int, String] = {
      def getPos(file: File): List[(Int, String)] = {
        val snpIterator = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(_.invalidSite).filter(_.isSNP)
        snpIterator.map(_ match {
          case SNP(r, c, a) => (c, r)
        }).toList
      }
      fList.flatMap(getPos(_)).toMap
    }

    /** Function to write in new text file. */
    def printToFile(f: File)(op: PrintWriter => Unit) {
      val writer = new PrintWriter(f)
      try { op(writer) } finally { writer.close() }
    }

    /** Return new name of samples with 9 characters and a space. */
    def truncateName(s: String): String = {
      if (s.length > 10) s.substring(s.length - 9, s.length) + " " //Cut s to length 10
      else { val res = "          ".substring(0, 10 - s.length); s + res } //Add spaces until length 10
    }

    /** List all VCFs from the path file, get all SNP positions, and print SNP sequence of each sample. */
    args.length match {
      case 2 => time {
        val fileList = Source.fromFile(new File(args(0))).getLines.map(new File(_)).toList
        val refMap = getPositions(fileList) // Map with ref. positions and bases
        println("Writing phy-file...")
        printToFile(new File(args(1))) { p =>
          p.println(fileList.size + 1 + " " + refMap.size) //Print total number of sequences (VCF's) + reference (1st sequence) and total number of SNP positions.
          p.print("H37RV_V5  ")
          refMap.keysIterator.toList.sorted.foreach(pos => p.print(refMap(pos)))
          fileList.foreach { file => // for each file print sequence
            val name = file.getParentFile.getName
            val snpMap = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(_.invalidSite).filter(_.isSNP).map(_ match {
              case SNP(r, c, a) => (c, a)
            }).toMap
            val nonSnpSet = (Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(_.invalidSite).filter(_.isSNP).filter(line =>
              refMap.contains(line.split("\t")(1).toInt)).map(line => line.split("\t")(1).toInt)).toSet
            val snpSeq = refMap.keysIterator.toList.sorted.map(pos =>
              if (snpMap.contains(pos)) snpMap(pos)
              else if (nonSnpSet.contains(pos)) "N"
              else refMap(pos)).mkString
            p.println
            p.print(truncateName(name) + snpSeq)
            println(name + ":\t" + snpMap.size + "\tSNPs")
          }
        }
        println("Total of " + fileList.size + " VCFs read.")
        println("Output: " + args(1))
      }
      case _ => println(usage)
    }
  }

}