package ara

import java.io.File
import scala.io.Source
import java.io.PrintWriter
import ara.Cluster._

/**
 * Extract SNPs from BlueJay output and convert them to 21 bp markers, given a reference genome.
 */
object BlueJayMarkers {
  
  case class Config(
    val bjOutput: File = null,
    val ref: File = null
  )
  
  def main(args: Array[String]) {
    
    val parser = new scopt.OptionParser[Config]("java -jar ara.jar bj-markers") {
      opt[File]("bj-output") required () action { (x, c) => c.copy(bjOutput = x) } text ("BlueJay output file to generate markers from lineage specific SNPs.")
      opt[File]("ref") required () action { (x, c) => c.copy(ref = x) } text { "Reference fasta file to extend SNPs with 10 bp." }
    }

    parser.parse(args, Config()) map { config =>
      
      /** Check if one of the two clusters contains the reference genome. */
      val clusters = config.bjOutput.getName.dropRight(16).split("_").toList
      println("Clusters: " + clusters.mkString(", "))
      
      println("One of the two clusters contains the reference genome: " + (clusters(0).hasReference || clusters(1).hasReference))      
      
      /** Reference genome to extend SNPs to 21 bp markers */
      val refGenome = Source.fromFile(config.ref).getLines.filterNot(_.startsWith(">")).mkString
      println("Reference genome size: " + refGenome.size)
      println
      
      /** Extract lineage associated SNPs from BlueJay output */
      val associatedSnps = Source.fromFile(config.bjOutput).getLines.filterNot(_.startsWith("#")).map{ line =>
        val arr = line.split("\t")
        (arr(10), arr(11), arr(12))      
      }.filter{m => 
        val arr = m._1.split("_")
        val ref = arr(1)
        val alt = arr(2)
        ref.size == 1 && alt.size == 1 // && m._2.equals("presence")
      }.toList.sortBy(_._1.split("_")(0).toInt)
      println(associatedSnps.size + " BlueJay SNPs")
      //associatedSnps.foreach(println)
      //println      
      
      /**
       * Remove SNP positions within 10 bp
       */
      val associatedSnpsPos = associatedSnps.map(_._1.split("_")(0).toInt).sorted //All SNP positions
      val associatedSnpsPos2 = associatedSnpsPos.filterNot { x =>
        val idx = associatedSnpsPos.indexOf(x)
        if (idx == 0) (associatedSnpsPos(idx + 1) - associatedSnpsPos(idx) < 11)
        else if (idx == associatedSnpsPos.size - 1) (associatedSnpsPos(idx) - associatedSnpsPos(idx - 1) < 11)
        else (associatedSnpsPos(idx) - associatedSnpsPos(idx - 1) < 11) || (associatedSnpsPos(idx + 1) - associatedSnpsPos(idx) < 11)
      }
      println(associatedSnpsPos2.size + " Non-overlapping SNPs")
      //associatedSnpsPos2.foreach(println)
      val associatedSnps2 = associatedSnps.filter(snp => associatedSnpsPos2.contains(snp._1.split("_")(0).toInt))
      //println(associatedSnps2.size)
      
      /** Make markers */
      def mkMarkers(snps: List[(String, String, String)]): List[(String, String)] = {
        snps.map{snp =>
        val arr = snp._1.split("_")
        val pos = arr(0).toInt
        val ref = arr(1)
        val alt = arr(2)
        if (ref != refGenome.substring(pos - 1, pos)) println("\tRef. base differs at position " + pos)
        val marker = refGenome.substring(pos - 11, pos - 1) + alt + refGenome.substring(pos, pos + 10)
        (">" + ref + pos + alt + "_" + snp._2 + "_" + snp._3, marker)
        }
      }      
      
      val markers = if (clusters(0).hasReference || clusters(1).hasReference) {
        mkMarkers(associatedSnps2.map{snp =>
          val cluster = snp._3
          val refCluster = clusters.filterNot(_.equals(cluster)).mkString
          val arr = snp._1.split("_")
          val coordinate = arr(0)
          val ref = arr(1)
          val alt = arr(2)
          val complementMarkers = List("A", "C", "G", "T").filterNot(_.equals(alt)).map(n => (coordinate + "_" + ref + "_" + n, "presence", refCluster))
          
          (snp) :: complementMarkers
        }.flatten)
      } else mkMarkers(associatedSnps2)
      println(markers.size + " SNP markers" + (if (clusters(0).hasReference || clusters(1).hasReference) " including complementary markers"))
      //markers.foreach(println)
      
      /**
       * Remove non-unique markers
       */
      val allMarkers = markers.map(_._2).toList      
      val mCounts = allMarkers.map(m1 => (m1, allMarkers.count(m2 => m2 == m1))).toMap
      //mCounts.foreach(println)
      val selection = markers.filter(snp => mCounts(snp._2) == 1)
      print(selection.size + " Unique markers")
      if (clusters(0).hasReference || clusters(1).hasReference) print(" including complementary markers")
      
      val pw = new PrintWriter(config.bjOutput.getName.dropRight(15) + "bluejaymarkers")
      pw.println("# BlueJay Lineage specific mutations")
      pw.println("# Console input: " + "scala BlueJayMarkers.scala " + args.mkString(" ") + "\n#")
      selection.foreach(m => pw.println(m._1 + "\n" + m._2))
      pw.close
    }
       
  }
}