//package macaw2

import scala.io.Source
import com.sun.org.apache.xalan.internal.xsltc.compiler.Sort
import java.io.PrintWriter
import java.io.File

/**
 * @author Arlin
 */
object Vcf2snpPhylip {
  val usage = "scala Vcf2snpPhylip.scala [directory] [output.phy] \nMultiple directories may be given as input and should be separated with spaces."

  def main(args: Array[String]) {

    /**
     * Elapsed time function
     */
    def time[R](block: => R): R = {
      val t0 = System.currentTimeMillis()
      val result = block // call-by-name
      val t1 = System.currentTimeMillis()
      println("Elapsed time: " + (t1 - t0) + "ms")
      result
    }

    /**
     * Function to read VCF-file and return a map with SNPs.
     * Key is position, value is a tuple (ref,alt).
     */
    def readVCF(f: File): Map[Int, (String, String)] = {
      val directories = f.toString().mkString.split("/") //For windows use ("""\\"""), for Linux ("/")
      val map1 = Map(0 -> (directories(directories.size - 2), directories(directories.size - 1))) //Store directory name and VCF name as tuple with key 0
      object SNP { //Object SNP to match with line in VCF.
        def unapply(s: String): Option[(Int, String, String)] = {
          val arr = s.mkString.split("\t")
          val r = arr(3)
          val a = arr(4)
          if (arr(6) == "PASS" && r.length() == 1 && a.length() == 1)
            Some((arr(1).toInt, r, a))
          else None
        }
      }
      def isSNP(line: String): Boolean = line match {
        case SNP(p, r, a) => true
        case _ => false
      }
      val lineIterator = Source.fromFile(f).getLines()
      val map2 = lineIterator.filterNot(_.startsWith("#")).filter(line => isSNP(line)).map(line => line match {
        case SNP(p, r, a) => (p -> (r, a))
      }).toMap
      map1 ++ map2
    }

    /**
     *  Write to file
     */
    def writeToFile(phyName: String, list: List[Map[Int, (String, String)]]) = {
      def printToFile(f: File)(op: PrintWriter => Unit) { //Function to write in new textfile. */
        val writer = new PrintWriter(f)
        try { op(writer) } finally { writer.close() }
      }
      def truncateName(s: (String, String)): String = {
        def truncate(name: String) = {
          if (name.length > 10) name.substring(name.length - 10, name.length - 1) //Cut name to length 10
          else { val res = "          ".substring(0, 10 - name.length); name + res } //Add spaces to length 10
        }
        if (s._2 == "reduced.vcf") truncate(s._1) //For Biek2012 and Blouin2012 files
        else truncate(s._2.substring(0, s._2.indexOf("."))) //Filename without extension
      }
      //Concatenate all VCF maps into a total reference map with all SNP positions and the reference base as value, and remove the value with key 0.
      val refMap: Map[Int, String] = list.flatMap(m => m.map(snp => (snp._1, snp._2._1))).toMap - 0
      printToFile(new File(phyName)) { p => //Use p to write in phy-file
        p.println(list.length + 1 + " " + refMap.size) //Print total number of sequences (VCF's) + reference (1st sequence) and total number of SNP positions.
        p.print("H37RV_V5  ")
        refMap.keysIterator.toList.sorted.foreach(pos => p.print(refMap(pos)))
        p.println
        for (snpMap <- list) {
          val nameVCF = truncateName(snpMap(0)) //Get name of the VCF
          println(nameVCF + ": " + snpMap.size + " SNPs") //Print in Console
          p.print(nameVCF) //Print in phy-file
          refMap.keysIterator.toList.sorted.foreach { pos =>
            if (!snpMap.contains(pos)) p.print(refMap(pos))
            else p.print(snpMap(pos)._2)
          }
          p.println
        }
      }
    }

    /**
     * List all VCF files in the given directory.
     */
    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles.toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.endsWith(".vcf")) => List(f)
      case _ => Nil
    }

    /**
     * Extract HQ SNPs, store in Map[Int, (String, String)] for each VCF and print information in a phylip format.
     */
    args.length match {
      case n if (n > 1) => time {
        val fileList = (0 to n - 2).toList.flatMap(idx => listFiles(new File(args(idx))))
        val mapList = fileList.map(file => readVCF(file))
        writeToFile(args(n - 1), mapList)
        println("Total of " + mapList.length + " files in all given directories.")
        println("Output: " + args(n - 1))
      }
      case _ => println(usage)
    }
  }

}