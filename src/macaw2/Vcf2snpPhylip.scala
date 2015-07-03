//package macaw2

import scala.io.Source
import com.sun.org.apache.xalan.internal.xsltc.compiler.Sort
import java.io.PrintWriter
import java.io.File

/**
 * This tool reads VCFs from the path file, and returns a SNP phy-file.
 */
object Vcf2snpPhylip {
  val usage = "scala Vcf2snpPhylip.scala [pathfile] [output.phy]"

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
    
    def invalidSite(line: String): Boolean = {
      val arr = line.split("\t")
      if (arr(4) == "." && arr(6) == "PASS") true
      else false
    }

    /** Get all SNP positions */
    def getPositions(fList: List[File]): Map[Int, String] = {
      def getPos(file: File): List[(Int, String)] = {
        val snpIterator = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(invalidSite(_)).filter(isSNP(_))
        snpIterator.map(_ match {
          case SNP(p, r, a) => (p, r)
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
        var excludesCount = 0
        var excludeList: List[(String, Int)] = List()
        println("Writing phy-file...")
        printToFile(new File(args(1))) { p =>
          p.println(fileList.size + 1 + " " + refMap.size) //Print total number of sequences (VCF's) + reference (1st sequence) and total number of SNP positions.
          p.print("H37RV_V5  ")
          refMap.keysIterator.toList.sorted.foreach(pos => p.print(refMap(pos)))
          fileList.foreach { file => // for each file print sequence
            val name = file.getParentFile.getName
            val snpMap = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filter(isSNP(_)).map( _ match {
              case SNP(p, r, a) => (p, a) 
            }).toMap
            val nonSnpSet = (Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(isSNP(_)).filter(line => 
              refMap.contains(line.split("\t")(1).toInt)).map(line => line.split("\t")(1).toInt)).toSet
            val snpSeq = refMap.keysIterator.toList.sorted.map(pos => 
              if (snpMap.contains(pos)) snpMap(pos)
              else if (nonSnpSet.contains(pos)) "N"
              else refMap(pos)
            ).mkString
            p.println
            p.print(truncateName(name) + snpSeq)
            println(name + ":\t" + snpMap.size + "\tSNPs")
            val ncount = snpSeq.filter(_ == 'N').size
            if (ncount.toFloat / refMap.size > 0.05) { //If number of N in sequence is more than 5% of total sequence.
              excludesCount = excludesCount + 1    
              excludeList = (name, ncount) :: excludeList
            }
          }
        }
        println("Total of " + fileList.size + " VCFs read.")
        println("WARNING: there are " + excludesCount + " sequences with more than 5% Ns.")
        println("Consider rerunning this tool without the samples:")
        println(excludeList.size)
        excludeList.foreach(x => println(x._1 + " has " + x._2 + "Ns."))
        println("Output: " + args(1))
      }
      case _ => println(usage)
    }
  }

}