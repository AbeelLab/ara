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
    val araResults: File = null,
    val trainingSet: File = null,
    val clusterDir: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar markers-qc") {
      opt[File]("markers") action { (x, c) => c.copy(markerFile = x) } text ("File containing marker sequences. This file has to be a multi-fasta file with the headers indicating the name of the markers.") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x + "_qc.txt") } text ("Output name.")
      opt[File]('d', "data-dir") required () action { (x, c) => c.copy(directory = x) } text ("Data directory.")
      opt[File]('r', "results") required () action { (x, c) => c.copy(araResults = x) } text ("Ara results file.")
      opt[File]('t', "trainingset") required () action { (x, c) => c.copy(trainingSet = x) } text ("Trainingset: VCF paths file.")
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
      /* Load SNP-typer files */
      val study = listAraFiles(config.directory)
      println("Ara SNP-typer outputs: " + study.size)
      /* Load Ara results */
      val predictions = tLines(config.araResults).map{line => val arr = line.split("\t"); (arr(1) -> arr(2))}.toMap
      val trnSet = tLines(config.trainingSet).toList.map(_.split("/")(8)).map(s => (s -> predictions(s))).toMap
      println("Training set size: " + trnSet.size)
      /* Load cluster files */
      val clusterFile = listClusterFiles(config.clusterDir).map(c => (c.getName.split("_")(0) -> c)).toMap
      val clusters = mtbcClusters.filterNot(_ == "MTBC").map(c => clusterFile(c)).sorted
      println("Clusters: " + clusters.size)
      println
      
      val pw = new PrintWriter(config.output)
      pw.println("$$\t" + (clusters.map(_.getName.dropRight(8)).sorted.mkString("\t")))
      markers.foreach { m =>
        
        pw.print(m)
        
        clusters.foreach{c => 
          
          val samples = tLines(c)
          
          pw.print(samples.size)
        }
        
        pw.println
      }
      pw.close

    }
  }
}