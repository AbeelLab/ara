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
    val bjOutput: File = null, val clusters: File = null,
    val ref: File = null, val filter: Boolean = false, val output: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar bj-markers") {
      opt[File]("bj-output") required () action { (x, c) => c.copy(bjOutput = x) } text ("BlueJay output directory or tsv-file to generate markers from lineage specific SNPs.")
      opt[File]("clusters") required () action { (x, c) => c.copy(clusters = x) } text ("Directory with cluster files.")
      opt[File]("ref") required () action { (x, c) => c.copy(ref = x) } text { "Reference fasta file to extend SNPs with 10 bp." }
      opt[Unit]("filter") action { (x, c) => c.copy(filter = true) } text ("Filter out groups of clusters with both clusters size 0, and descendants.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Prefix output file.")
    }

    parser.parse(args, Config()) map { config =>

      val tsvFiles = config.bjOutput.listFiles().filter(_.getName.endsWith("_lineage_snp.tsv")).map(tsv => (tsv.getName.dropRight(16), tsv)).toMap
      val clusterFiles = config.clusters.listFiles().filter(_.getName.endsWith("_cluster")).map(c => (c.getName.dropRight(8), c)).toMap
      val tsvNames = tsvFiles.keySet
      val clusterNames = clusterFiles.keySet
      val fileSet = tsvNames.intersect(clusterNames).map(set => (set, tsvFiles(set))).toList.sortBy(_._1)

      def getClusterSize(cluster: String): Int = {
        Source.fromFile(clusterFiles(cluster)).getLines.filterNot(_.startsWith("#")).map { line => val arr = line.split("\t"); (arr(0), arr(1)) }.filterNot(_._2.startsWith("not")).size
      }

      def readTsvFile(file: File): List[(String, String, String)] = {
        Source.fromFile(file).getLines.filterNot(_.startsWith("#")).map { line =>
          val arr = line.split("\t")
          (arr(10), arr(11), arr(12))
        }.filter { m =>
          val arr = m._1.split("_")
          val ref = arr(1)
          val alt = arr(2)
          ref.size == 1 && alt.size == 1 && !m._3.startsWith("not")
        }.toList.sortBy(_._1.split("_")(0).toInt) // Sort SNPs by coordinate
      }

      def getSnpNumbers(set: (String, File)): (String, Map[String, Any]) = set match {
        case (name, tsv) => {
          /** Group by presence/absence and map size of associated SNP-list */
          val snps = readTsvFile(tsv).groupBy(_._2).mapValues(_.size)
          val presence = if (snps.contains("presence")) snps("presence") else 0
          val absence = if (snps.contains("absence")) snps("absence") else 0
          (name, Map("presenceSNPs" -> presence, "absenceSNPs" -> absence))
        }
        case _ => (null, Map.empty[String, Any])
      }

      /** Get number of SNPs and sample size, and group by ancestor. */
      val numbers = fileSet.map(getSnpNumbers(_)).groupBy(_._1.getAncestor).toList
      println(numbers.flatMap(_._2).size + " clusters")
      /** Filter out clusters that do not have a secondary group. */
      val numbers2 = numbers.filter(_._2.size == 2)
      println(numbers2.flatMap(_._2).size + " Grouped clusters")

      /** List clusters that should be left out, because this cluster and their companion cluster both have 0 associated SNPs. */
      val emptyGroups = numbers2.filter {
        _ match {
          case (ancestor, children) => {
            val c1Numbers = children(0)._2
            val c2Numbers = children(1)._2
            c1Numbers("presenceSNPs") == 0 && c1Numbers("absenceSNPs") == 0 && c2Numbers("presenceSNPs") == 0 && c2Numbers("absenceSNPs") == 0
          }
        }
      }.flatMap(_._2).map(_._1)

      /** Remove clusters listed in emptyGroups and their descendants. */
      println("Clusters left out: ")
      val filteredClusters = numbers2.flatMap(_._2).filter {
        _ match {
          case (cluster, map) => {
            val matches = emptyGroups.filter(c => cluster.contains(c)) // Cluster matches with emptyGroups
            if (!matches.isEmpty) println("\t" + cluster)
            matches.isEmpty
          }
        }
      }.map(_._1)
      println(filteredClusters.size + " clusters in total")

      /** Reference genome to extend SNPs to 21 bp markers */
      val refGenome = Source.fromFile(config.ref).getLines.filterNot(_.startsWith(">")).mkString
      println("Reference genome " + config.ref.getName + " size: " + refGenome.size)
      println

      /** Extract lineage associated SNPs from BlueJay output for valid clusters */
      val associatedSnps = filteredClusters.map(c => (c, readTsvFile(tsvFiles(c)))).toMap // Map[cluster, snps]      
      val totalSnps = associatedSnps.flatMap(_._2).size
      println(totalSnps + " cluster-specific SNPs from all clusters.")
      //associatedSnps.foreach(println)
      //val L1to4Snps = associatedSnps.flatMap(_._2).filter(_._3 == "L1-L2-L3-L4").size
      //println(L1to4Snps + " SNPs from L1-L2-L3-L4")
      //println(totalSnps - L1to4Snps + " SNPs, L1-L2-L3-L4 excluded")
      //val associatedSnps2 = associatedSnps.filterNot(_._1 == "L1-L2-L3-L4")
      //val overlap = associatedSnps2.flatMap(_._2).groupBy(_._1).filter(_._2.size > 1) //Group by mutation and count occurences
      //println(totalSnps - L1to4Snps - overlap.size + " Unique SNPs")

      def removeWithin10bp(ls: List[(String, String, String)]): List[(String, String, String)] = {
        val allPositions = ls.map(_._1.split("_")(0).toInt).toList.sorted
        val allPositions2 = if (allPositions.size < 2) allPositions else allPositions.filterNot { x =>
          val idx = allPositions.indexOf(x)
          if (idx == 0) (allPositions(idx + 1) - allPositions(idx) < 11)
          else if (idx == allPositions.size - 1) (allPositions(idx) - allPositions(idx - 1) < 11)
          else (allPositions(idx) - allPositions(idx - 1) < 11) || (allPositions(idx + 1) - allPositions(idx) < 11)
        }
        ls.filter(snp => allPositions2.contains(snp._1.split("_")(0).toInt))
      }

      /** Remove SNP positions within 10 bp of another SNP within the cluster */
      //val emptyL1toL4 = Map(("L1-L2-L3-L4", List.empty[(String, String, String)]))
      val associatedSnps3 = associatedSnps.map(_ match {
        case (cluster, snps) => (cluster, removeWithin10bp(snps))
      })// ++ emptyL1toL4
      println(associatedSnps3.flatMap(_._2).size + " total SNPs not within 10bp of another SNP in the same cluster")

      /** Make markers */
      def mkMarkers(snps: List[(String, String, String)]): List[(String, String)] = {
        snps.map { snp =>
          val arr = snp._1.split("_")
          val pos = arr(0).toInt
          val ref = arr(1)
          val alt = arr(2)
          if (ref != refGenome.substring(pos - 1, pos)) println("\tRef. base differs at position " + pos)
          val marker = refGenome.substring(pos - 11, pos - 1) + alt + refGenome.substring(pos, pos + 10)
          (">" + ref + pos + alt + "_" + snp._2 + "_" + snp._3, marker)
        }
      }

      def mkComplementMarkers(snps: List[(String, String, String)]): List[(String, String)] = {
        val complements = snps.map { snp =>
          val mutation = snp._1.split("_")
          val coordinate = mutation(0)
          val ref = mutation(1)
          val alt = mutation(2)
          val complementMarkers = List("A", "C", "G", "T").filterNot(_ == alt).map(n => (coordinate + "_" + ref + "_" + n, "presence", snp._3))
          complementMarkers
        }.flatten
        mkMarkers(complements)
      }

      /** Markers for each cluster grouped by ancestor*/
      /**val markers = associatedSnps3.groupBy(_._1.getAncestor).toList.map {
        _ match {
          case (ancestor, clusters) => {
            val c1 = clusters.toList(0)
            val c2 = clusters.toList(1)
            val complementMarkers = List("A", "C", "G", "T")
            val mLists = if (c1._2.isEmpty) {
              Map((c1._1, mkComplementMarkers(c1._1, c2._2)), (c2._1, mkMarkers(c2._2)))
            } else if (c2._2.isEmpty) {
              Map((c1._1, mkMarkers(c1._2)), (c2._1, mkComplementMarkers(c2._1, c1._2)))
            } else {
              clusters.map(c => (c._1, mkMarkers(c._2)))
            }
            mLists
          }
        }
      }.flatten.toMap*/
      val markers = associatedSnps3.map{_ match {
        case (clusters, snps) => {
          if (snps.map(_._2).contains("absence")) (clusters, mkComplementMarkers(snps))
          else (clusters, mkMarkers(snps))
        }
      }}
      println(markers.flatMap(_._2).size + " markers, incl. complementary markers")

      def filterMarkers(ls: List[(String, String)]): List[(String, String)] = {
        val markers = ls.map(_._2)
        val mCounts = markers.map(m1 => (m1, markers.count(m2 => m2 == m1))).toMap
        ls.filter(m => mCounts(m._2) == 1)
      } 
      
      /** Remove non-unique markers */
      val selection = markers.map {
        _ match {
          case (cluster, markers) => (cluster, filterMarkers(markers))
        }
      }
      println(selection.flatMap(_._2).size + " Unique markers")

      val mpw = new PrintWriter(config.output + "_bluejaymarkers")
      mpw.println("# BlueJay Lineage specific mutations")
      mpw.println("# Console input: " + "scala BlueJayMarkers.scala " + args.mkString(" ") + "\n#")
      selection.flatMap(_._2).toList.sortBy(_._1.split("_")(0).drop(2).dropRight(1).toInt).foreach(m => mpw.println(m._1 + "\n" + m._2))
      mpw.close

      /** Print numbers */

      def printAllNumbers(ls: List[String]): Unit = {
        val groupedClusters = ls.groupBy(_.getAncestor).toList.sortBy(_._1)
        val npw = new PrintWriter(config.output + "_numbers.txt")
        npw.println("#Cluster\tSample size\tPresence SNPs\tAbsence SNPs\tFiltered presence SNPs\tFiltered absence SNPs\tPresence markers\tAbsence markers\tFiltered presence markers\tFiltered absence markers")

        groupedClusters.foreach {
          _ match {
            case (ancestor, clusters) => clusters.foreach { c =>
              val sampleSize = getClusterSize(c)
              val presenceSNPs = associatedSnps(c).filter(_._2 == "presence").size
              val absenceSNPs = associatedSnps(c).filter(_._2 == "absence").size
              val filteredPresenceSNPs = associatedSnps3(c).filter(_._2 == "presence").size
              val filteredAbsenceSNPs = associatedSnps3(c).filter(_._2 == "absence").size
              val presenceMarkers = markers(c).filter(_._1.split("_")(1) == "presence").size
              val absenceMarkers = markers(c).filter(_._1.split("_")(1) == "absence").size
              val filteredPresenceMarkers = selection(c).filter(_._1.split("_")(1) == "presence").size
              val filteredAbsenceMarkers = selection(c).filter(_._1.split("_")(1) == "absence").size

              npw.println(c + "\t" + sampleSize + "\t" + presenceSNPs + "\t" + absenceSNPs + "\t" + filteredPresenceSNPs + "\t" + filteredAbsenceSNPs + "\t" +
                presenceMarkers + "\t" + absenceMarkers + "\t" + filteredPresenceMarkers + "\t" + filteredAbsenceMarkers) // + "\t" + numbers("clusterPath") + numbers("tsvPath"))
            }
          }
        }
        npw.close
      }

      printAllNumbers(filteredClusters)

    }

  }
}