package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object CountBJsnps {

  case class Config(val bjOutputs: File = null, val clusters: File = null, val markers: File = null, filter: Boolean = false, val out: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar count-bj-snps") {
      opt[File]("bj-outputs") required () action { (x, c) => c.copy(bjOutputs = x) } text ("Directory with BlueJay output files with cluster specific SNPs.")
      opt[File]("clusters") required () action { (x, c) => c.copy(clusters = x) } text ("Directory with cluster files.")
      opt[File]("markers") required () action { (x, c) => c.copy(markers = x) } text ("Directory with marker files.")
      opt[Unit]("filter") action { (x, c) => c.copy(filter = true) } text ("Filter out groups of clusters with both clusters size 0, and descendants.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(out = x + ".txt") } text ("Prefix output file.")
    }

    parser.parse(args, Config()) map { config =>

      val tsvFiles = config.bjOutputs.listFiles().filter(_.getName.endsWith("_lineage_snp.tsv")).map(tsv => (tsv.getName.dropRight(16), tsv)).toMap
      val clusterFiles = config.clusters.listFiles().filter(_.getName.endsWith("_clusters")).map(c => (c.getName.dropRight(9), c)).toMap
      val markerFiles = config.markers.listFiles().filter(_.getName.endsWith("_bluejaymarkers")).map(m => (m.getName.dropRight(15), m)).toMap

      val tsvNames = tsvFiles.keySet
      val clusterNames = clusterFiles.keySet
      val markerNames = markerFiles.keySet

      val fileSet = tsvNames.intersect(clusterNames).intersect(markerNames).map(set => (set, tsvFiles(set), clusterFiles(set), markerFiles(set))).toList.sortBy(_._1)

      def getNumbers(set: (String, File, File, File)): (String, List[(String, Map[String, Any])]) = set match {
        case (name, tsv, cluster, markers) => {

          /** Group by cluster, group by presence/absence and map size of associated SNP-list */
          val snps = Source.fromFile(tsv).getLines.filterNot(_.startsWith("#")).map { line =>
            val arr = line.split("\t")
            (arr(10), arr(11), arr(12))
          }.toList.groupBy(_._3).map(c => (c._1, c._2.groupBy(_._2).mapValues(_.size)))

          /** Sample size of cluster */
          val size = Source.fromFile(cluster).getLines.filterNot(_.startsWith("#")).map { line => val arr = line.split("\t"); (arr(1), arr(0)) }.toList.groupBy(_._1).mapValues(_.size)

          /** Group by cluster, group by presence/absence and map size of marker list */
          val markerSet = Source.fromFile(markers).getLines.filterNot(_.startsWith("#")).grouped(2).map(m => m(0)).toList.groupBy(_.split("_")(2)).map(c => (c._1, c._2.groupBy(_.split("_")(1)).mapValues(_.size)))

          val map = name.split("_").map { c =>

            val presenceSNPs = if (snps.contains(c)) { if (snps(c).contains("presence")) snps(c)("presence") else 0 } else 0
            val absenceSNPs = if (snps.contains(c)) { if (snps(c).contains("absence")) snps(c)("absence") else 0 } else 0
            val presenceMarkers = if (markerSet.contains(c)) { if (markerSet(c).contains("presence")) markerSet(c)("presence") else 0 } else 0
            val absenceMarkers = if (markerSet.contains(c)) { if (markerSet(c).contains("absence")) markerSet(c)("absence") else 0 } else 0
            (c, Map("cluster" -> c, "sampleSize" -> size(c), "presenceMarkers" -> presenceMarkers, "absenceMarkers" -> absenceMarkers, "presenceSNPs" -> presenceSNPs, "absenceSNPs" -> absenceSNPs, "clusterPath" -> cluster.getCanonicalPath, "markerPath" -> markers.getCanonicalPath, "tsvPath" -> tsv.getCanonicalPath))
          }.toList
          (name, map)
        }
        case _ => (null, List.empty[(String, Map[String, Any])])
      }

      val pw = new PrintWriter(config.out)
      pw.println("#Number of presence and absence SNPs in BlueJay tsv-file output.")
      pw.println("#Cluster\tSample size\tPresence markers\tAbsence markers\tPresence SNPs\tAbsence SNPs\tCluster file-path\tMarker file-path\tTsv file-path")
      val numbers = fileSet.map(getNumbers(_))
      numbers.foreach { cg =>
        cg match {
          case (cGroup, clusters) => clusters.foreach {
            _ match {
              case (c, numbers) => pw.println(c + "\t" + numbers("sampleSize") + "\t" + numbers("presenceMarkers") + "\t" + numbers("absenceMarkers") + "\t" + numbers("presenceSNPs") + "\t" + numbers("absenceSNPs") + "\t" + numbers("clusterPath") + "\t" + numbers("markerPath") + "\t" + numbers("tsvPath"))
            }
          }
        }
      }
      pw.close

      if (config.filter) {

        val emptyGroups = numbers.filter {
          _ match {
            case (cGroup, clusters) => {
              val c1Numbers = clusters(0)._2
              val c2Numbers = clusters(1)._2
              c1Numbers("presenceMarkers") == 0 && c1Numbers("absenceMarkers") == 0 && c2Numbers("presenceMarkers") == 0 && c2Numbers("absenceMarkers") == 0
            }
          }
        }.flatMap(_._2).map(_._1)
        
        println("Groups of clusters left out: ")
        val filteredGroups = numbers.filter{
          _ match {
            case (cGroup, clusters) => {
              val matches = emptyGroups.filter(c => cGroup.contains(c))
              if (!matches.isEmpty) println("\t" + cGroup)
              matches.isEmpty
            }
          }

        }
        
        
        val pw2 = new PrintWriter("allMarkers.txt")
        pw2.println("# All markers for each hierarchical clusters.")
        pw2.println("# Command: java -jar ara.jar all-bj-markers " + args.mkString(" "))
        val filteredFiles = filteredGroups.map(g => markerFiles(g._1))
        val allMarkers = filteredFiles.flatMap(file => Source.fromFile(file).getLines.filterNot(_.startsWith("#")).grouped(2))
        pw2.println("# Cluster groups: " + filteredFiles.size)
        pw2.println("# Total markers: " + allMarkers.size)
        allMarkers.foreach(m => pw2.println(m.mkString("\n")))
        pw2.close

      }

    }

  }
}