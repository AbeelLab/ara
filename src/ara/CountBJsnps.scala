package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object CountBJsnps {

  case class Config(val bjOutputs: File = null, val clusters: File = null, val out: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar count-bj-snps") {
      opt[File]("bj-outputs") required () action { (x, c) => c.copy(bjOutputs = x) } text ("Directory with BlueJay output files with cluster specific SNPs.")
      opt[File]("clusters") required () action { (x, c) => c.copy(clusters = x) } text ("Directory with cluster files.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(out = x + ".txt") } text ("Prefix output file.")
    }

    parser.parse(args, Config()) map { config =>

      val tsvFiles = config.bjOutputs.listFiles().filter(_.getName.endsWith("_lineage_snp.tsv")).map(tsv => (tsv.getName.dropRight(16), tsv)).toMap
      val clusterFiles = config.clusters.listFiles().filter(_.getName.endsWith("_clusters")).map(c => (c.getName.dropRight(9), c)).toMap
      
      val tsvNames = tsvFiles.keySet
      val clusterNames = clusterFiles.keySet
      val fileSet = tsvNames.intersect(clusterNames).map(set => (set, tsvFiles(set), clusterFiles(set)))
      
      
      
      val pw = new PrintWriter(config.out)
      pw.println("#Number of presence and absence SNPs in BlueJay tsv-file output.")
      pw.println("#Cluster\tSample size\tPresence\tAbsence\tTsv file-path\tCluster file-path")
      fileSet.foreach { set => set match {
        case (name, tsv, cluster) => {
            val snps = Source.fromFile(tsv).getLines.filterNot(_.startsWith("#")).map { line =>
              val arr = line.split("\t")
              (arr(10), arr(11), arr(12))
            }.toList.groupBy(_._3).map(c => (c._1, c._2.groupBy(_._2).mapValues(_.size)))// Group by cluster, group by presence/absence and map size of markerlist
            val size = Source.fromFile(cluster).getLines.filterNot(_.startsWith("#")).map {line => val arr = line.split("\t"); (arr(1), arr(0))}.toList.groupBy(_._1).mapValues(_.size)
            
            name.split("_").foreach{ c =>              
              if (snps.contains(c)) { // If cluster has markers
                val markers = snps(c)
                pw.println(c + "\t" + size(c) + "\t" + (if (markers.contains("presence")) markers("presence") else 0) + "\t" + (if (markers.contains("absence")) markers("absence") else 0) + "\t" + tsv + "\t" + cluster)
              }
              else pw.println(c + "\t" + size(c) + "\t" + 0 + "\t" + 0 + "\t" + tsv + "\t" + cluster)
            }
            
        }
      }
        
              
      }

      pw.close

    }

  }
}