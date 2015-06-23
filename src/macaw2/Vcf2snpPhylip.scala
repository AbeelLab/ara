package macaw2

import scala.io.Source
import com.sun.org.apache.xalan.internal.xsltc.compiler.Sort
import java.io.PrintWriter
import java.io.File

/**
 * This tool only reads VCFs that are named reduced.vcf
 */
object Vcf2snpPhylip {
  val usage = "scala Vcf2snpPhylip.scala [directory] [output.phy] \nMultiple directories may be given as input and should be separated with spaces."

  def main(args: Array[String]) {

    /** Elapsed time function */
    def time[R](block: => R): R = {
      val t0 = System.currentTimeMillis()
      val result = block // call-by-name
      val t1 = System.currentTimeMillis()
      println("Elapsed time: " + (t1 - t0) + "ms")
      result
    }

    /** Object SNP to match with line in VCF. */
    object SNP {
      def unapply(s: String): Option[(Int, String, String)] = {
        val arr = s.mkString.split("\t")
        val r = arr(3)
        val a = arr(4)
        if (arr(6) == "PASS" && r.length() == 1 && a.length() == 1)
          Some((arr(1).toInt, r, a))
        else None
      }
    }

    /** Determine if line represents a SNP. */
    def isSNP(line: String): Boolean = line match {
      case SNP(p, r, a) => true
      case _ => false
    }

    /** Get all SNP positions */
    def getPositions(fList: List[File]): Set[(Int, String)] = {
      def getPos(file: File): List[(Int, String)] = {
        val snpIterator = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filter(isSNP(_))
        snpIterator.map(_ match {
          case SNP(p, r, a) => (p, r)
        }).toList
      }
      fList.flatMap(getPos(_)).toSet
    }

    /** Function to write in new textfile. */
    def printToFile(f: File)(op: PrintWriter => Unit) {
      val writer = new PrintWriter(f)
      try { op(writer) } finally { writer.close() }
    }
    
    /** Return new name of samples with 9 characters. */
    def truncateName(s: String): String = {
      if (s.length > 10) s.substring(s.length - 9, s.length) + " " //Cut s to length 10
      else { val res = "          ".substring(0, 10 - s.length); s + res } //Add spaces until length 10
    }

     /** List all VCF files in the given directory. */
    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles.toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.equals("reduced.vcf")) => List(f)
      case _ => Nil
    }

    /** List all VCFs from the given directories, get all SNP positions, and print SNP sequence of each sample. */
    args.length match {
      case n if (n > 1) => time {
        val fileList = (0 to n - 2).toList.flatMap(idx => listFiles(new File(args(idx))))
        val refList = getPositions(fileList).toList.sorted
        printToFile(new File(args(n - 1))) { p =>
          p.println(fileList.size + 1 + " " + refList.size) //Print total number of sequences (VCF's) + reference (1st sequence) and total number of SNP positions.
          p.print("H37RV_V5  ")
          refList.foreach(i => p.print(i._2))
          p.println
          fileList.foreach { file =>
            val name = file.getParentFile.getName
            p.print(truncateName(name))
            val snpMap = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filter(isSNP(_)).map( _ match {
              case SNP(p, r, a) => (p, a) 
            }).toMap
            val nonSnpSet = (Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(isSNP(_)).map(line => line.split("\t")(1).toInt)).toSet
            println(name + ":\t" + snpMap.size + "\tSNPs,\t" + nonSnpSet.size + "\tnon-SNP variants")
            val snpSeq = (refList.map(pos =>
              if (snpMap.contains(pos._1)) snpMap(pos._1)
              else if (nonSnpSet.contains(pos._1)) "N"
              else pos._2)).mkString
            p.println(snpSeq)
          }
        }
        println("Total of " + fileList.size + " files in all given directories.")
        println("Output: " + args(n - 1))
      }
      case _ => println(usage)
    }
  }

}