package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import ara.Cluster._
import atk.util.Tool

object QCGlobalMarkers extends MTBCclusters with Tool {

  case class Config(
    val markerFile: File = null,
    val output: String = null,
    val directory: File = null,
    val clusterDir: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar markers-qc") {
      opt[File]("markers") action { (x, c) => c.copy(markerFile = x) } text ("File containing marker sequences. This file has to be a multi-fasta file with the headers indicating the name of the markers.") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x + "_qc.txt") } text ("Output name.")
      opt[File]('d', "data-dir") required () action { (x, c) => c.copy(directory = x) } text ("Data directory.")
      opt[File]('c', "cluster-dir") required () action { (x, c) => c.copy(clusterDir = x) } text ("Directory with cluster files.")
    }

    parser.parse(args, Config()) map { config =>

      def listAraFiles(f: Any): List[File] = f match {
        case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listAraFiles(_))
        case f: File if (f.isFile() && f.getName.endsWith(".ara") && !f.getName.endsWith("interpret.ara")) => List(f)
        case _ => Nil
      }

      def listClusterFiles(f: Any): List[File] = f match {
        case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listClusterFiles(_))
        case f: File if (f.isFile() && f.getName.endsWith("_cluster")) => List(f)
        case _ => Nil
      }

      /* Load markers */
      val lines = if (config.markerFile != null) tLines(config.markerFile).toList else scala.io.Source.fromInputStream(MacawSNPtyper.getClass().getResourceAsStream("/hierclusters_global_bluejaymarkers")).getLines().filterNot(f => f.startsWith("#") || f.trim.size == 0).toList
      val markers = lines.grouped(2).map(f => (f(0).substring(1))).toList.sortBy(_.split("_")(2))
      println("Markers: " + markers.size)
      /* Load SNP-typer files */
      val araFile = listAraFiles(config.directory).map(a => (a.getName.dropRight(4) -> a)).toMap
      println("Ara SNP-typer outputs: " + araFile.size)
      /* Load cluster files */
      val clusterFile = listClusterFiles(config.clusterDir).map(c => (c.getName.split("_")(0) -> c)).toMap
      val clusters = markers.map(_.split("_")(2)).distinct.sorted
      println("Clusters with  markersets: " + clusters.size)
      println

      def markerDetected(marker: String, sample: String): Boolean = {
        if (araFile.keySet.contains(sample)) {
          val markerInFile = tLines(araFile(sample)).dropRight(2).map {
            _ match {
              case Marker(m, c, p) => new ClusterMarker(m, c, p)
            }
          }.filter(m => m.mutInfo == marker).head
          //println(marker + " is present in " + sample + ": " + markerInFile.isPresent)
          markerInFile.isPresent
        } else {
          false
        }        
      }

      val pw = new PrintWriter(config.output)
      pw.println("$$\t" + (clusters.mkString("\t")))
      markers.foreach { m =>

        pw.print(m + "\t")
        println(m)
        
        clusters.map(c => clusterFile(c)).foreach { c =>
          
          val samples = tLines(c).filterNot(_.contains("not")).map(_.split("\t")(0))
          val present =  samples.map(s => markerDetected(m, s)).count(_ == true)
          //print("\t" + c  + "\t" + present)
          //println
          pw.print(present + "\t")
        }

        pw.println
      }
      pw.close

    }
  }
}