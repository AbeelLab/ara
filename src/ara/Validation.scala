package ara

import java.io.PrintWriter
import scala.io.Source
import be.abeel.bioinformatics.FastaIterator
import atk.util.Tool
import java.io.File
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileWriterFactory
import net.sf.samtools.SAMRecordIterator
import edu.northwestern.at.utils.math.randomnumbers.MersenneTwister

object Validation extends Tool {

  override val description = """This program randomly subsets two bam files to have have  a simulated mixed infection"""

  case class Config(
    val input: Seq[File] = Seq(),
    val output: String = null,
    val count: Int = 10000000,
    val probs: Seq[Int] = Seq(),
    val replicates: Int = 1)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar generate-validation") {
      opt[String]("prefix") required () action { (x, c) => c.copy(output = x) } text ("Output prefix")
      opt[Seq[File]]('f', "files") required () action { (x, c) => c.copy(input = x) } text ("List of input bam files")
      opt[Int]("count") action { (x, c) => c.copy(count = x) } text ("Approx. total number of reads in output. Default=10,000,000 ")
      opt[Seq[Int]]('p', "percentages") required () action { (x, c) => c.copy(probs = x) } text ("Percentages of mixes.")
      opt[Int]('r', "replicates") action { (x, c) => c.copy(replicates = x)} text ("Number of replicates. Default=1")
    }

    parser.parse(args, Config()) map { config =>

      assume(config.input.size == 2, "This program requires two input bam files")

      val sam1 = new SAMFileReader(config.input(0));
      val sam2 = new SAMFileReader(config.input(1));

      val seqs1 = sam1.getFileHeader().getSequenceDictionary().size()
      println("references: " + seqs1)

      val seqs2 = sam2.getFileHeader().getSequenceDictionary().size()
      println("references: " + seqs2)
      assume(sam1.hasIndex(), "The bam files needs to be indexed things...")
      assume(sam2.hasIndex(), "The bam files needs to be indexed things...")

      val totalReads1 =
        (for (i <- 0 until seqs1)
          yield sam1.getIndex().getMetaData(i).getAlignedRecordCount).sum
      val totalReads2 =
        (for (i <- 0 until seqs2)
          yield sam2.getIndex().getMetaData(i).getAlignedRecordCount).sum

      println("File 1: " + config.input(0) + "\t" + totalReads1)
      println("File 2: " + config.input(1) + "\t" + totalReads2)

      val probSeq = config.probs

      println("Probability table:")
      val probTable = for (i <- probSeq) yield {
        val prob1 = (i / 100.0) * (config.count.toDouble / totalReads1)
        val estimated1 = prob1 * totalReads1
        val prob2 = ((100 - i) / 100.0) * (config.count.toDouble / totalReads2)
        val estimated2 = prob2 * totalReads2
        println(i + "\t" + (100 - i) + "\t" + prob1 + "\t" + prob2 + "\t" + estimated1.toInt + "\t" + estimated2.toInt + "\t" + (estimated1 + estimated2).toInt)

        ((i, 100 - i, prob1, prob2))

      }

      def createMixes(n: Int) = {

        for (r <- 1 to n) {
          println("-----replicate " + r + "-----")
          val outputs = for (i <- probSeq) yield {

            val outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(sam1.getFileHeader(),
              false, new File(config.output + "." + "r" + r + "." + i + ".bam"));

            outputSam
          }

          val workMatrix = probTable.zip(outputs)

          /* Process first file*/
          var counter = 0
          val rg = new MersenneTwister
          println("Processing 1")
          val it: SAMRecordIterator = sam1.iterator()

          while (it.hasNext()) {
            val record = it.next()
            counter += 1
            if (counter % 100000 == 0)
              print(".")
            for ((prob, out) <- workMatrix) {
              val r = rg.nextDouble()
              if (r <= prob._3)
                out.addAlignment(record)
            }
          }
          it.close()
          println
          /* Process second file */
          val it2: SAMRecordIterator = sam2.iterator()
          println("Processing 2")
          while (it2.hasNext()) {
            val record = it2.next()
            counter += 1
            if (counter % 100000 == 0)
              print(".")
            for ((prob, out) <- workMatrix) {
              val r = rg.nextDouble()
              if (r <= prob._4)
                out.addAlignment(record)
            }
          }
          it2.close()
          for (oo <- outputs) {
            println("Closing: " + oo)
            oo.close
          }
        }

      }

      createMixes(config.replicates)

      println("Finished.")
    }

  }
}