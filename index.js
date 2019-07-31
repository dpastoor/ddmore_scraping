const Nightmare = require('nightmare')

const http = require('http');
const fs = require('fs');
const  _ = require('lodash');
const leftPad = require('left-pad');

const BASELINE_URL = "http://repository.ddmore.eu/model";
const START_INDEX = 234;
const LAST_INDEX = 301;

let download = function(url, dest, cb) {
  var file = fs.createWriteStream(dest);
  var request = http.get(url, function(response) {
    response.pipe(file);
    file.on('finish', function() {
      file.close(cb);  // close() is async, call cb after close completes.
    });
  }).on('error', function(err) { // Handle errors
    fs.unlink(dest); // Delete the file async. (But we don't check the result)
    if (cb) cb(err.message);
  });
};

let getFileNamesAndVersion = function() {
    // Scrape the links from the leaf dropdown showing what files available
    var links = document.querySelectorAll('li.jstree-leaf a');
    // title attribute corresponds to file name
    let files = Array.prototype.map.call(links, function (e) {
        return e.getAttribute('title')
    });
    // multiple bulleted versions, of which the final version is appended
    // to the download path, so want to either grab the latest version or default
    // to 1
    let version = document.querySelector('#History ul:last-child').innerText
    .split('\n')
    .reduce((acc, line) => {
      if (line.includes("Version") & acc == 1) {
          return parseInt(line.split(":")[1].trim()) // line looks like: "Version: <num>"
      }
      return acc
    }, 1)

    return {files, version};
  }

let downloadFiles = function(i) {
  const paddedVidIndex = leftPad(i, 3, 0);
  console.log(`${BASELINE_URL}/DDMODEL00000${paddedVidIndex}#Files`)
  new Nightmare()
    .goto(`${BASELINE_URL}/DDMODEL00000${paddedVidIndex}#Files`)
    .wait()
    .click('#treeView ul li:nth-child(1) ul li a span')
    .wait('#filegoeshere')
    .evaluate(getFileNamesAndVersion)
    .run((err, res) =>{
      console.log(res)
      console.log('about to download');
      if (!res) {
        fs.appendFile('missingModels.txt', paddedVidIndex +"\n", (err) => console.log(err))
      } else {
        _.forEach(res.files, (file) => {
          console.log('downloading: ' + file)
          download(`${BASELINE_URL}/download/DDMODEL00000${paddedVidIndex}.${res.version}?filename=${file}`, paddedVidIndex + "/" + file, console.log)
        })
          fs.writeFile(paddedVidIndex + "/" + paddedVidIndex + '.json', JSON.stringify(res, null, 4), (err) => console.log(err))
          console.log('done downloading');
      }
    })
}

let i = START_INDEX;
let dlInterval = setInterval(() => {
  if (i > LAST_INDEX) {
    clearInterval(dlInterval)
    return
  }
  let dir = leftPad(i, 3, 0);
  if (!fs.existsSync(dir)) {
    fs.mkdirSync(dir)
  }
  downloadFiles(i)
  i++
}, 5000)

