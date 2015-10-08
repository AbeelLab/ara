package ara

import java.io.File
import scala.io.Source
import java.io.PrintWriter

/**
 * Extract SNPs from BlueJay output and convert them to 21 bp markers, given a reference genome.
 */
object BlueJayMarkers {
  
  val usage = "scala BlueJayMarkers.scala [BlueJay output] [reference fasta] [output]"
  
  def main(args: Array[String]) {
    if (args.size != 3) println(usage)
    else {
      
      val associatedSnps = Source.fromFile(new File(args(0))).getLines.filterNot(_.startsWith("#")).map{ line =>
        val arr = line.split("\t")
        (arr(10), arr(11), arr(12))      
      }.filter{m => 
        val arr = m._1.split("_")
        val ref = arr(1)
        val alt = arr(2)
        ref.size == 1 && alt.size == 1
      }.toList
      
      val refGenome = Source.fromFile(new File(args(1))).getLines.filterNot(_.startsWith(">")).mkString
      println("refGenome.size: " + refGenome.size)
      
      val markers = associatedSnps.map{snp =>
        val arr = snp._1.split("_")
        val pos = arr(0).toInt
        val ref = arr(1)
        val alt = arr(2)
        println("Ref. base matches at position " + pos + ": " + (ref == refGenome.substring(pos - 1, pos)))
        val marker = refGenome.substring(pos - 11, pos - 1) + alt + refGenome.substring(pos, pos + 10)
        (">" + ref + pos + alt + "_" + snp._2 + "_" + snp._3, marker)
      }
      
      val pw = new PrintWriter(new File(args(2)))
      pw.println("# BlueJay Lineage specific mutations")
      pw.println("# Console input: " + "scala BlueJayMarkers.scala " + args.mkString(" ") + "\n#")
      markers.foreach(m => pw.println(m._1 + "\n" + m._2))
      pw.close
    }
    
  }
}